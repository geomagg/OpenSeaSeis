/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/*
 * Module INTERFEROMETRIA: virtual shot gather (VSG) from ambient noise by cross-correlation with a
 * reference trace (virtual source), stacked over time windows and over the whole input.
 * Author: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026.
 * Same algorithm as the SeaView plugin "Interferometria" (cseis.plugin.csPluginInterferometria).
 *
 * Input: ensembles = one time chunk with all nodes (e.g. 5-min traces after CONCATENATE, sorted by
 *        time_samp1 and with ENS_DEFINE header time_samp1). Without ensembles, the whole input is one chunk.
 * Per chunk, per trace and per window of 'window' seconds:
 *   demean, 5% Tukey taper, band-pass (cosine-tapered mask fmin-fmax),
 *   temporal normalisation (none, one-bit, running absolute mean), spectral whitening (optional),
 *   then conj(R) X is accumulated in the frequency domain, R = reference trace with the same 'pair' header.
 * At the end of the input: C(tau) = IFFT( sum conj(R) X ), tau = -max_lag..max_lag (lag 0 at sample nl),
 * or C(tau)+C(-tau), tau = 0..max_lag ('sides sum'). One output trace per node ('node' header) and pair value.
 * tau > 0: the trace lags the reference.
 */

#include "cseis_includes.h"
#include "csStandardHeaders.h"
#include <cmath>
#include <cstring>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

namespace mod_interferometria {
  struct Acc {
    double pairVal, nodeVal;
    vector<double> re, im;     // sum over windows of conj(R) X  [nfft]
    int numWin;
    csTrace* trace;            // output trace (header copied from the first input trace of this node)
    double x, y;               // node coordinates
  };
  struct VariableStruct {
    int hdrId_ref, hdrId_pair, hdrId_node;
    int hdrId_t1, hdrId_t1us, hdrId_x, hdrId_y;
    int hdrId_off, hdrId_fold, hdrId_lag0, hdrId_sides;
    double refValue;
    int nsIn; double dt;
    int nw, nl, hop, nwin, nfft, nout, nram;
    double fmin, fmax;
    int tempNorm;              // 0 none, 1 one-bit, 2 running absolute mean
    bool whiten, symmetric, normalize, setOffset;
    vector<double>* mask;
    vector<double>* taper;
    vector<Acc*>* accs;
    map< pair<double,double>, int >* accIndex;
    map< double, pair<double,double> >* refXY;
    csTraceGather* outGather;
    int numEnsembles, numNoRef, numTimeMismatch, numDeadWin, numNoRefTraces;
  };
  static void fft( double* re, double* im, int n, bool inverse );
  static void preprocess( float const* s, int i0, VariableStruct const* v, double* re, double* im );
}
using namespace mod_interferometria;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_interferometria_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csTraceHeaderDef* hdef = env->headerDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars   = new VariableStruct();
  edef->setVariables( vars );
  edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  vars->mask = NULL; vars->taper = NULL; vars->accs = NULL; vars->accIndex = NULL; vars->refXY = NULL; vars->outGather = NULL;
  vars->numEnsembles = vars->numNoRef = vars->numTimeMismatch = vars->numDeadWin = vars->numNoRefTraces = 0;
  string text;

  //--- Virtual source
  param->getString( "ref", &text, 0 );
  if( !hdef->headerExists( text ) ) writer->error( "Header '%s' (parameter 'ref') does not exist", text.c_str() );
  vars->hdrId_ref = hdef->headerIndex( text );
  param->getDouble( "ref", &vars->refValue, 1 );
  string refName = text;

  //--- Pairing header (each trace is correlated with the reference that has the same value, e.g. chan)
  vars->hdrId_pair = -1;
  text = hdef->headerExists( "chan" ) ? "chan" : "none";
  if( param->exists( "pair" ) ) param->getString( "pair", &text );
  if( text.compare( "none" ) ) {
    if( !hdef->headerExists( text ) ) writer->error( "Header '%s' (parameter 'pair') does not exist", text.c_str() );
    vars->hdrId_pair = hdef->headerIndex( text );
  }
  string pairName = text;

  //--- Node identification (one output trace per node and pair value)
  text = hdef->headerExists( "rcv" ) ? "rcv" : refName;
  if( param->exists( "node" ) ) param->getString( "node", &text );
  if( !hdef->headerExists( text ) ) writer->error( "Header '%s' (parameter 'node') does not exist", text.c_str() );
  vars->hdrId_node = hdef->headerIndex( text );
  string nodeName = text;

  //--- Times
  vars->nsIn = shdr->numSamples;
  vars->dt   = shdr->sampleInt / 1000.0;
  double maxLag = 10.0, winLen = 0.0, overlap = 0.0;
  if( param->exists( "max_lag" ) ) param->getDouble( "max_lag", &maxLag );
  if( param->exists( "window" ) )  param->getDouble( "window", &winLen );
  if( param->exists( "overlap" ) ) param->getDouble( "overlap", &overlap );
  if( maxLag <= 0 ) writer->error( "max_lag must be > 0" );
  if( overlap < 0 || overlap > 90 ) writer->error( "overlap must be between 0 and 90 %%" );
  int ns = vars->nsIn;
  int nw = ( winLen > 0 ) ? std::min( ns, (int)lround( winLen / vars->dt ) ) : ns;
  nw = std::max( nw, 4 );
  int nl = std::min( (int)lround( maxLag / vars->dt ), nw - 1 );
  int hop = std::max( 1, (int)lround( nw * ( 1.0 - overlap / 100.0 ) ) );
  if( nw >= ns ) hop = nw;
  int nwin = std::max( 1, 1 + ( ns - nw ) / hop );
  int nfft = 1;
  while( nfft < nw + nl ) nfft *= 2;
  vars->nw = nw; vars->nl = nl; vars->hop = hop; vars->nwin = nwin; vars->nfft = nfft;

  //--- Band, normalisation
  double fnyq = 0.5 / vars->dt;
  vars->fmin = 0.0; vars->fmax = fnyq;
  if( param->exists( "band" ) ) {
    param->getDouble( "band", &vars->fmin, 0 );
    if( param->getNumValues( "band" ) > 1 ) param->getDouble( "band", &vars->fmax, 1 );
  }
  if( vars->fmax <= 0 || vars->fmax > fnyq ) vars->fmax = fnyq;
  if( vars->fmax <= vars->fmin ) writer->error( "band: fmax must be > fmin" );
  vars->tempNorm = 2;
  double ramWin = 30.0;
  if( param->exists( "temp_norm" ) ) {
    param->getString( "temp_norm", &text );
    if( !text.compare( "none" ) ) vars->tempNorm = 0;
    else if( !text.compare( "onebit" ) ) vars->tempNorm = 1;
    else if( !text.compare( "ram" ) ) vars->tempNorm = 2;
    else writer->error( "Unknown option for 'temp_norm': '%s'", text.c_str() );
    if( param->getNumValues( "temp_norm" ) > 1 ) param->getDouble( "temp_norm", &ramWin, 1 );
  }
  if( ramWin <= 0 ) writer->error( "temp_norm: running window must be > 0" );
  vars->nram = std::max( 1, (int)lround( ramWin / vars->dt ) );
  vars->whiten = true;
  if( param->exists( "whitening" ) ) { param->getString( "whitening", &text ); vars->whiten = !text.compare( "yes" ); }
  vars->symmetric = false;
  if( param->exists( "sides" ) ) {
    param->getString( "sides", &text );
    if( !text.compare( "sum" ) ) vars->symmetric = true;
    else if( text.compare( "two" ) ) writer->error( "Unknown option for 'sides': '%s'", text.c_str() );
  }
  vars->normalize = true;
  if( param->exists( "normalize" ) ) { param->getString( "normalize", &text ); vars->normalize = !text.compare( "yes" ); }

  //--- Headers
  vars->hdrId_t1   = hdef->headerExists( HDR_TIME_SAMP1.name ) ? hdef->headerIndex( HDR_TIME_SAMP1.name ) : -1;
  vars->hdrId_t1us = hdef->headerExists( HDR_TIME_SAMP1_US.name ) ? hdef->headerIndex( HDR_TIME_SAMP1_US.name ) : -1;
  vars->hdrId_x = hdef->headerExists( "rec_x" ) ? hdef->headerIndex( "rec_x" ) : -1;
  vars->hdrId_y = hdef->headerExists( "rec_y" ) ? hdef->headerIndex( "rec_y" ) : -1;
  vars->setOffset = ( vars->hdrId_x >= 0 && vars->hdrId_y >= 0 );
  if( param->exists( "offset" ) ) {
    param->getString( "offset", &text );
    if( !text.compare( "no" ) ) vars->setOffset = false;
    else if( vars->setOffset == false ) writer->error( "offset yes: headers rec_x and rec_y are required" );
  }
  vars->hdrId_off = -1;
  if( vars->setOffset ) {
    if( !hdef->headerExists( HDR_OFFSET.name ) ) hdef->addStandardHeader( HDR_OFFSET.name );
    vars->hdrId_off = hdef->headerIndex( HDR_OFFSET.name );
  }
  if( !hdef->headerExists( HDR_FOLD.name ) ) hdef->addStandardHeader( HDR_FOLD.name );
  vars->hdrId_fold = hdef->headerIndex( HDR_FOLD.name );
  if( !hdef->headerExists( "lag0_ms" ) ) hdef->addHeader( TYPE_FLOAT, "lag0_ms", "Time of lag 0 in the correlation trace [ms] (INTERFEROMETRIA)" );
  vars->hdrId_lag0 = hdef->headerIndex( "lag0_ms" );
  if( !hdef->headerExists( "xc_sides" ) ) hdef->addHeader( TYPE_INT, "xc_sides", "2: lags -max..max, 1: both sides summed, lags 0..max (INTERFEROMETRIA)" );
  vars->hdrId_sides = hdef->headerIndex( "xc_sides" );

  //--- Masks
  int nfft2 = nfft;
  double df = 1.0 / ( nfft2 * vars->dt );
  vars->mask = new vector<double>( nfft2, 0.0 );
  double w = std::max( 2*df, 0.1 * ( vars->fmax - vars->fmin ) );
  for( int f = 0; f <= nfft2/2; f++ ) {
    double fr = f * df, v;
    if( fr < vars->fmin - w || fr > vars->fmax + w ) v = 0.0;
    else if( fr < vars->fmin ) v = ( vars->fmin <= 0 ) ? 1.0 : 0.5 - 0.5 * cos( M_PI * ( fr - ( vars->fmin - w ) ) / w );
    else if( fr > vars->fmax ) v = 0.5 + 0.5 * cos( M_PI * ( fr - vars->fmax ) / w );
    else v = 1.0;
    if( f == 0 ) v = 0.0;    // no DC
    (*vars->mask)[f] = v;
    if( f > 0 && f < nfft2/2 ) (*vars->mask)[nfft2 - f] = v;
  }
  vars->taper = new vector<double>( nw, 1.0 );
  int ne = std::max( 1, (int)( 0.05 * nw ) );
  for( int i = 0; i < nw; i++ ) {
    if( i < ne ) (*vars->taper)[i] = 0.5 - 0.5 * cos( M_PI * i / ne );
    else if( i >= nw - ne ) (*vars->taper)[i] = 0.5 - 0.5 * cos( M_PI * ( nw - 1 - i ) / ne );
  }
  vars->accs = new vector<Acc*>();
  vars->accIndex = new map< pair<double,double>, int >();
  vars->refXY = new map< double, pair<double,double> >();
  vars->outGather = new csTraceGather( hdef );

  vars->nout = vars->symmetric ? nl + 1 : 2 * nl + 1;
  shdr->numSamples = vars->nout;

  static char const* NORMS[3] = { "none", "one-bit", "running absolute mean" };
  writer->line( "  Virtual source: %s = %g   pairing header: %s   node header: %s", refName.c_str(), vars->refValue, pairName.c_str(), nodeName.c_str() );
  char buf[128];
  buf[0] = '\0';
  if( vars->tempNorm == 2 ) snprintf( buf, sizeof(buf), " (%g s)", ramWin );
  writer->line( "  Band: %g - %g Hz   temporal normalisation: %s%s   whitening: %s",
                vars->fmin, vars->fmax, NORMS[vars->tempNorm], buf, vars->whiten ? "yes" : "no" );
  writer->line( "  Window: %d samples (%g s), %d window(s) per input trace, start every %d samples (%g s); FFT %d",
                nw, nw * vars->dt, nwin, hop, hop * vars->dt, nfft );
  buf[0] = '\0';
  if( !vars->symmetric ) snprintf( buf, sizeof(buf), ", lag 0 at %g ms", nl * vars->dt * 1000.0 );
  writer->line( "  Max lag: %g s (%d samples)   output: %s, %d samples%s",
                nl * vars->dt, nl, vars->symmetric ? "both sides summed, lags 0..max" : "lags -max..max", vars->nout, buf );
  writer->line( "  Offset header = distance to the virtual source: %s", vars->setOffset ? "yes" : "no (rec_x/rec_y missing or offset no)" );
  int unused = ns - ( ( nwin - 1 ) * hop + nw );
  if( unused * 10 > ns ) writer->warning( "The last %g s of each input trace are not used (window/overlap do not divide the trace)", unused * vars->dt );
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_interferometria_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csTraceHeaderDef const* hdef = env->headerDef;
  int ntr = traceGather->numTraces();
  int nfft = vars->nfft;

  //---------------------------------------------------------------------------
  // One time chunk: accumulate conj(R) X for all traces
  if( ntr > 0 ) {
    vars->numEnsembles++;
    // References: every trace of the virtual source in this chunk (there may be several files/times)
    vector<int> refs;
    vector<double> tStart( ntr, 0.0 );
    for( int i = 0; i < ntr; i++ ) {
      csTraceHeader const* h = traceGather->trace(i)->getTraceHeader();
      if( vars->hdrId_t1 >= 0 ) tStart[i] = h->intValue( vars->hdrId_t1 ) + ( vars->hdrId_t1us >= 0 ? h->intValue( vars->hdrId_t1us ) * 1.0e-6 : 0.0 );
      if( fabs( h->doubleValue( vars->hdrId_ref ) - vars->refValue ) > 1.0e-6 * std::max( 1.0, fabs( vars->refValue ) ) ) continue;
      refs.push_back( i );
      double pv = ( vars->hdrId_pair >= 0 ) ? h->doubleValue( vars->hdrId_pair ) : 0.0;
      if( vars->setOffset && vars->refXY->find( pv ) == vars->refXY->end() ) {
        (*vars->refXY)[pv] = make_pair( h->doubleValue( vars->hdrId_x ), h->doubleValue( vars->hdrId_y ) );
      }
    }
    map<int,int> refOf;       // used references (trace index -> trace index), for the window loop
    if( refs.empty() ) {
      vars->numNoRef++;
      if( vars->numNoRef <= 5 ) writer->warning( "Time chunk %d (%d traces): no trace of the virtual source, chunk skipped", vars->numEnsembles, ntr );
    }
    else {
      // Accumulator and reference of each trace: same pairing value AND same start time (within dt/2)
      vector<int> accOf( ntr, -1 ), refTrc( ntr, -1 );
      for( int i = 0; i < ntr; i++ ) {
        csTraceHeader const* h = traceGather->trace(i)->getTraceHeader();
        double pv = ( vars->hdrId_pair >= 0 ) ? h->doubleValue( vars->hdrId_pair ) : 0.0;
        int rk = -1;
        bool samePair = false;
        double bestDt = 1.0e30;
        for( size_t k = 0; k < refs.size(); k++ ) {
          csTraceHeader const* hr = traceGather->trace( refs[k] )->getTraceHeader();
          double pr = ( vars->hdrId_pair >= 0 ) ? hr->doubleValue( vars->hdrId_pair ) : 0.0;
          if( pr != pv ) continue;
          samePair = true;
          double d = fabs( tStart[i] - tStart[ refs[k] ] );
          if( d < bestDt ) { bestDt = d; rk = refs[k]; }
        }
        if( !samePair ) { vars->numNoRefTraces++; continue; }
        if( bestDt > 0.5 * vars->dt ) {
          vars->numTimeMismatch++;
          if( vars->numTimeMismatch <= 5 ) writer->warning( "Chunk %d: trace %d: no reference trace starting at the same time (closest %.4f s): not correlated. Use ALIGN_TIME", vars->numEnsembles, i+1, bestDt );
          continue;
        }
        refOf[rk] = rk;
        double nv = h->doubleValue( vars->hdrId_node );
        pair<double,double> key( pv, nv );
        map< pair<double,double>, int >::iterator it = vars->accIndex->find( key );
        int ia;
        if( it == vars->accIndex->end() ) {
          Acc* a = new Acc();
          a->pairVal = pv; a->nodeVal = nv;
          a->re.assign( nfft, 0.0 ); a->im.assign( nfft, 0.0 );
          a->numWin = 0;
          a->trace = vars->outGather->createTrace( hdef, vars->nout );
          a->trace->getTraceHeader()->copyFrom( h );
          a->x = ( vars->hdrId_x >= 0 ) ? h->doubleValue( vars->hdrId_x ) : 0.0;
          a->y = ( vars->hdrId_y >= 0 ) ? h->doubleValue( vars->hdrId_y ) : 0.0;
          ia = (int)vars->accs->size();
          vars->accs->push_back( a );
          (*vars->accIndex)[key] = ia;
        }
        else ia = it->second;
        accOf[i] = ia;
        refTrc[i] = rk;
      }
      // Window by window: reference spectra first, then every trace
      vector<double> xr( nfft ), xi( nfft );
      map< int, vector<double> > rRe, rIm;
      map< int, bool > rDead;
      for( int w = 0; w < vars->nwin; w++ ) {
        int i0 = w * vars->hop;
        for( map<int,int>::const_iterator r = refOf.begin(); r != refOf.end(); ++r ) {
          vector<double>& re = rRe[r->second];
          vector<double>& im = rIm[r->second];
          re.resize( nfft ); im.resize( nfft );
          float const* s = traceGather->trace( r->second )->getTraceSamples();
          bool dead = true;
          for( int k = 0; k < vars->nw; k++ ) if( s[i0+k] != 0.0f ) { dead = false; break; }
          rDead[r->second] = dead;
          if( !dead ) preprocess( s, i0, vars, &re[0], &im[0] );
        }
        for( int i = 0; i < ntr; i++ ) {
          if( accOf[i] < 0 ) continue;
          int k = refTrc[i];
          if( rDead[k] ) { vars->numDeadWin++; continue; }
          float const* s = traceGather->trace(i)->getTraceSamples();
          bool dead = true;
          for( int j = 0; j < vars->nw; j++ ) if( s[i0+j] != 0.0f ) { dead = false; break; }
          if( dead ) { vars->numDeadWin++; continue; }
          if( i == k ) { xr = rRe[k]; xi = rIm[k]; }
          else preprocess( s, i0, vars, &xr[0], &xi[0] );
          Acc* a = (*vars->accs)[ accOf[i] ];
          double const* rr = &rRe[k][0];
          double const* ri = &rIm[k][0];
          for( int f = 0; f < nfft; f++ ) {            // conj(R) * X
            a->re[f] += rr[f] * xr[f] + ri[f] * xi[f];
            a->im[f] += rr[f] * xi[f] - ri[f] * xr[f];
          }
          a->numWin++;
        }
      }
    }
    traceGather->freeAllTraces();
  }

  if( !edef->isLastCall() ) {
    edef->setTracesAreWaiting();    // output comes at the end of the input
    return;
  }

  //---------------------------------------------------------------------------
  // End of input: inverse FFT and output
  int nacc = (int)vars->accs->size();
  if( nacc == 0 ) {
    writer->warning( "No correlation computed (virtual source not found?)" );
    return;
  }
  int nl = vars->nl;
  vector<double> re( nfft ), im( nfft ), c( 2*nl+1 );
  for( int ia = 0; ia < nacc; ia++ ) {
    Acc* a = (*vars->accs)[ia];
    re = a->re; im = a->im;
    fft( &re[0], &im[0], nfft, true );
    for( int j = -nl; j <= nl; j++ ) c[j+nl] = re[ ( j + nfft ) % nfft ];
    float* o = a->trace->getTraceSamples();
    if( vars->symmetric ) for( int j = 0; j <= nl; j++ ) o[j] = (float)( c[nl+j] + c[nl-j] );
    else for( int j = 0; j < 2*nl+1; j++ ) o[j] = (float)c[j];
    if( vars->normalize ) {
      float m = 0.0f;
      for( int j = 0; j < vars->nout; j++ ) m = std::max( m, fabsf( o[j] ) );
      if( m > 0.0f ) for( int j = 0; j < vars->nout; j++ ) o[j] /= m;
    }
    csTraceHeader* h = a->trace->getTraceHeader();
    h->setIntValue( vars->hdrId_fold, a->numWin );
    h->setFloatValue( vars->hdrId_lag0, vars->symmetric ? 0.0f : (float)( nl * vars->dt * 1000.0 ) );
    h->setIntValue( vars->hdrId_sides, vars->symmetric ? 1 : 2 );
    if( vars->setOffset ) {
      map< double, pair<double,double> >::const_iterator r = vars->refXY->find( a->pairVal );
      if( r != vars->refXY->end() ) h->setFloatValue( vars->hdrId_off, (float)hypot( a->x - r->second.first, a->y - r->second.second ) );
    }
    a->re.clear(); a->im.clear();
  }
  vars->outGather->moveTracesTo( 0, vars->outGather->numTraces(), traceGather );
  writer->line( "INTERFEROMETRIA: %d time chunk(s), %d output correlation(s)", vars->numEnsembles, nacc );
  if( vars->numNoRef > 0 )        writer->line( "  %d chunk(s) without the virtual source (skipped)", vars->numNoRef );
  if( vars->numNoRefTraces > 0 )  writer->line( "  %d trace(s) without a reference with the same pairing value", vars->numNoRefTraces );
  if( vars->numTimeMismatch > 0 ) writer->line( "  %d trace(s) not correlated: start time differs from the reference", vars->numTimeMismatch );
  if( vars->numDeadWin > 0 )      writer->line( "  %d window(s) skipped: zero samples (gaps) in the trace or in the reference", vars->numDeadWin );
  for( int ia = 0; ia < nacc && ia < 2000; ia++ ) {
    Acc* a = (*vars->accs)[ia];
    if( a->numWin == 0 ) writer->line( "  node %g (pair %g): no window stacked, output is zero", a->nodeVal, a->pairVal );
  }
}

//*************************************************************************************************
// Pre-processing of one window [i0, i0+nw) -> spectrum (length nfft)
//*************************************************************************************************
namespace mod_interferometria {
static void preprocess( float const* s, int i0, VariableStruct const* v, double* re, double* im ) {
  int nfft = v->nfft, nw = v->nw;
  vector<double> const& taper = *v->taper;
  vector<double> const& mask  = *v->mask;
  double mean = 0.0;
  for( int i = 0; i < nw; i++ ) mean += s[i0 + i];
  mean /= nw;
  for( int i = 0; i < nfft; i++ ) { re[i] = ( i < nw ) ? ( s[i0 + i] - mean ) * taper[i] : 0.0; im[i] = 0.0; }
  fft( re, im, nfft, false );
  for( int f = 0; f < nfft; f++ ) { re[f] *= mask[f]; im[f] *= mask[f]; }
  if( v->tempNorm != 0 ) {
    fft( re, im, nfft, true );
    if( v->tempNorm == 1 ) {
      for( int i = 0; i < nfft; i++ ) re[i] = ( i < nw ) ? ( re[i] > 0 ? 1.0 : ( re[i] < 0 ? -1.0 : 0.0 ) ) * taper[i] : 0.0;
    }
    else {
      vector<double> a( nw + 1, 0.0 ), o( nw );
      for( int i = 0; i < nw; i++ ) a[i + 1] = a[i] + fabs( re[i] );
      int h = v->nram / 2;
      for( int i = 0; i < nw; i++ ) {
        int j0 = std::max( 0, i - h ), j1 = std::min( nw, i + h + 1 );
        double m = ( a[j1] - a[j0] ) / ( j1 - j0 );
        o[i] = ( m > 0 ) ? re[i] / m : 0.0;
      }
      for( int i = 0; i < nfft; i++ ) re[i] = ( i < nw ) ? o[i] * taper[i] : 0.0;
    }
    for( int i = 0; i < nfft; i++ ) im[i] = 0.0;
    fft( re, im, nfft, false );
    for( int f = 0; f < nfft; f++ ) { re[f] *= mask[f]; im[f] *= mask[f]; }
  }
  if( v->whiten ) {
    // divide by the amplitude spectrum smoothed over a few bins
    int h = nfft / 2 + 1;
    int half = std::max( 1, nfft / 400 );
    vector<double> c( h + 1, 0.0 ), sm( h );
    for( int f = 0; f < h; f++ ) c[f + 1] = c[f] + hypot( re[f], im[f] );
    double big = 0.0;
    for( int f = 0; f < h; f++ ) {
      int f0 = std::max( 0, f - half ), f1 = std::min( h, f + half + 1 );
      sm[f] = ( c[f1] - c[f0] ) / ( f1 - f0 );
      big = std::max( big, sm[f] );
    }
    double eps = 1.0e-6 * big;
    for( int f = 0; f < h; f++ ) {
      double g = mask[f] / ( sm[f] + eps );
      re[f] *= g; im[f] *= g;
      if( f > 0 && f < nfft / 2 ) { re[nfft - f] *= g; im[nfft - f] *= g; }
    }
  }
}

/** In-place radix-2 complex FFT, exp(-i...) forward; inverse with 1/n */
static void fft( double* re, double* im, int n, bool inverse ) {
  for( int i = 1, j = 0; i < n; i++ ) {
    int bit = n >> 1;
    for( ; j & bit; bit >>= 1 ) j ^= bit;
    j ^= bit;
    if( i < j ) { std::swap( re[i], re[j] ); std::swap( im[i], im[j] ); }
  }
  for( int len = 2; len <= n; len <<= 1 ) {
    double ang = 2.0 * M_PI / len * ( inverse ? 1.0 : -1.0 );
    double wr = cos( ang ), wi = sin( ang );
    int half = len >> 1;
    for( int i = 0; i < n; i += len ) {
      double cr = 1.0, ci = 0.0;
      for( int k = 0; k < half; k++ ) {
        int a = i + k, b = a + half;
        double tr = re[b] * cr - im[b] * ci, ti = re[b] * ci + im[b] * cr;
        re[b] = re[a] - tr; im[b] = im[a] - ti;
        re[a] += tr;        im[a] += ti;
        double t = cr * wr - ci * wi;
        ci = cr * wi + ci * wr; cr = t;
      }
    }
  }
  if( inverse ) for( int i = 0; i < n; i++ ) { re[i] /= n; im[i] /= n; }
}
}

//*************************************************************************************************
// Parameter definition
//*************************************************************************************************
void params_mod_interferometria_( csParamDef* pdef ) {
  pdef->setModule( "INTERFEROMETRIA", "Virtual shot gather (VSG) from ambient noise",
    "Cross-correlation of every trace with the virtual source (reference trace with the same 'pair' header value), "
    "stacked over time windows and over the whole input. Each ensemble must be one time chunk with all nodes "
    "(e.g. ENS_DEFINE header fileno after INPUT_HDF5); without ensembles the whole input is one chunk. "
    "An ensemble may also hold several time chunks: each trace is correlated with the virtual-source trace that starts at the same time (within dt/2). "
    "Output: one correlation per node and pair value, at the end of the input. Times in seconds. "
    "Same algorithm as the SeaView plugin 'Interferometria'." );

  pdef->addParam( "ref", "Virtual source: header name and value", NUM_VALUES_FIXED );
  pdef->addValue( "rcv", VALTYPE_STRING, "Header name (e.g. rcv, node)" );
  pdef->addValue( "", VALTYPE_NUMBER, "Header value of the virtual source" );

  pdef->addParam( "pair", "Pairing header", NUM_VALUES_FIXED, "Each trace is correlated with the reference trace that has the same value of this header (e.g. chan = component). 'none': one reference for all traces" );
  pdef->addValue( "chan", VALTYPE_STRING, "Header name, or 'none'" );

  pdef->addParam( "node", "Node header", NUM_VALUES_FIXED, "Traces with the same node and pair values are stacked into one output trace" );
  pdef->addValue( "rcv", VALTYPE_STRING, "Header name" );

  pdef->addParam( "band", "Frequency band [Hz]", NUM_VALUES_VARIABLE, "Cosine-tapered mask (taper width 10% of the band)" );
  pdef->addValue( "0", VALTYPE_NUMBER, "fmin [Hz]" );
  pdef->addValue( "", VALTYPE_NUMBER, "fmax [Hz] (default: Nyquist)" );

  pdef->addParam( "temp_norm", "Temporal normalisation", NUM_VALUES_VARIABLE );
  pdef->addValue( "ram", VALTYPE_OPTION );
  pdef->addOption( "none", "No temporal normalisation" );
  pdef->addOption( "onebit", "One-bit normalisation (sign of the band-passed signal)" );
  pdef->addOption( "ram", "Running absolute mean normalisation" );
  pdef->addValue( "30", VALTYPE_NUMBER, "Running window [s] (option ram), e.g. half of the longest period of interest or longer" );

  pdef->addParam( "whitening", "Spectral whitening", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Divide each window spectrum by its smoothed amplitude spectrum" );
  pdef->addOption( "no", "No whitening" );

  pdef->addParam( "max_lag", "Maximum correlation lag [s]", NUM_VALUES_FIXED );
  pdef->addValue( "10", VALTYPE_NUMBER, "Maximum lag [s]" );

  pdef->addParam( "window", "Stacking window [s]", NUM_VALUES_FIXED, "Each input trace is cut into windows of this length; correlations of all windows and chunks are stacked" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Window length [s]. 0 = whole input trace" );

  pdef->addParam( "overlap", "Overlap between consecutive windows [%]", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "0 to 90" );

  pdef->addParam( "sides", "Output sides", NUM_VALUES_FIXED );
  pdef->addValue( "two", VALTYPE_OPTION );
  pdef->addOption( "two", "Lags -max_lag..max_lag; lag 0 at sample max_lag (header lag0_ms)" );
  pdef->addOption( "sum", "Sum of both sides C(t)+C(-t), lags 0..max_lag" );

  pdef->addParam( "normalize", "Normalise each output correlation (max |C| = 1)", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Normalise" );
  pdef->addOption( "no", "Keep the stacked amplitude" );

  pdef->addParam( "offset", "Set header offset = distance to the virtual source", NUM_VALUES_FIXED, "Requires rec_x and rec_y" );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Set offset" );
  pdef->addOption( "no", "Do not change offset" );
}

//************************************************************************************************
// Start exec phase
//*************************************************************************************************
bool start_exec_mod_interferometria_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

//************************************************************************************************
// Cleanup phase
//*************************************************************************************************
void cleanup_mod_interferometria_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( vars->accs != NULL ) {
    for( size_t i = 0; i < vars->accs->size(); i++ ) delete (*vars->accs)[i];
    delete vars->accs;
  }
  if( vars->accIndex != NULL ) delete vars->accIndex;
  if( vars->refXY != NULL ) delete vars->refXY;
  if( vars->outGather != NULL ) delete vars->outGather;
  if( vars->mask != NULL ) delete vars->mask;
  if( vars->taper != NULL ) delete vars->taper;
  delete vars;
}

extern "C" void _params_mod_interferometria_( csParamDef* pdef ) {
  params_mod_interferometria_( pdef );
}
extern "C" void _init_mod_interferometria_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_interferometria_( param, env, writer );
}
extern "C" bool _start_exec_mod_interferometria_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_interferometria_( env, writer );
}
extern "C" void _exec_mod_interferometria_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_interferometria_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_interferometria_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_interferometria_( env, writer );
}
