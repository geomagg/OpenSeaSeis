/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module SPK_LSQMATCH: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include "csSeismicReader.h"
#include <cstring>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_LSQMATCH
 *
 * SEISPAK process lsqmatch (xlsqmatch.f / elsqmatch.f), trace-wise types 1, 3, 5:
 * least-squares (Wiener) match filter of ONLINE2 to ONLINE1, computed and applied
 * trace by trace, on the whole trace (TINWL = 0, LSQSUB0) or in sliding time
 * windows (TINWL, TOUTWL, TOUTOV, LSQSUB) with trapezoid merging (MTRAPZ).
 *   filter: autocorrelation of ONLINE2 (LEN lags, zero lag x (1+PDIAG/100)) and
 *   cross-correlation ONLINE1 x ONLINE2 starting at lag PRED (PRED = -LEN/2:
 *   two-sided match filter), solved by Levinson (LEVI).
 *   TYPE 1: out = ONLINE1 - F*ONLINE2   TYPE 3: out = ONLINE1 + F*ONLINE2
 *   TYPE 5: out = F*ONLINE2
 * ONLINE1 is the flow; ONLINE2 is read from the file given in 'filename' (.cseis),
 * trace by trace in the same order. Optional moveout shift before the design
 * (T0, VEL: t = sqrt(T0^2 + (1000*RANG/VEL)^2)) and back after it. If an ONLINE1
 * trace is zero (max < .001) the ONLINE2 trace is zeroed, as in the original.
 * The filter is written to FILTER (Fortran unformatted, dvsout format read by
 * mfilt2 / SPK_MFILT2): record 1 = NF NTR NREC SR(=1.), record 2 = NF coefficients.
 *
 * SEASEIS notes:
 *  - xlsqmatch wrote every filter (one per trace and window) and then rewrote FILTER
 *    with the FIRST one only; filter first (default) does the same, filter mean
 *    writes the average of all filters.
 *  - Types 2, 4, 6 (multi-trace, record-wise) and the X windows are not ported yet.
 *  - The header of FILTER has NTR = number of filters computed, NREC = 1.
 */

namespace mod_spk_lsqmatch {
  struct VariableStruct {
    int   type;            // 1, 3, 5
    int   nflt, iprd;      // filter length, prediction distance (samples)
    float pdiag;
    int   liw, low, lov;   // windows in samples (elm 4,5,6)
    double t0, vel;
    int   idOff;
    int   nt; float sr;
    csSeismicReader* reader;
    csTraceHeaderDef* hdef2;
    char* hdrBlock;
    int   nt2;
    std::vector<float> z;
    std::string filterFile;
    bool  filterMean;
    std::vector<double> fsum; std::vector<float> ffirst; int nfilt;
    int   ntrZeroZ, ntrIer;
  };
}
using namespace mod_spk_lsqmatch;

//---------------------------------------------------------------- numeric routines (1-based, as Fortran)
// CORRELG: C(LA:LB) = A correl B ; c points to C(LA)
static void correlg( float const* A, int N, float const* B, int M, float* c, int LA, int LB ) {
  for( int i = LA; i <= LB; i++ ) c[i-LA] = 0.0f;
  for( int I = std::max(1, LA-N+1); I <= std::min(M, LB); I++ ) {
    int MK = std::max(I, LA), K = MK-I+1, IB = M-I+1;
    for( int J = MK; J <= std::min(N+I-1, LB); J++ ) { c[J-LA] = B[IB-1]*A[K-1] + c[J-LA]; K++; }
  }
}
// CONVG: C(LA:LB) = A conv B
static void convg( float const* A, int N, float const* B, int M, float* c, int LA, int LB ) {
  for( int i = LA; i <= LB; i++ ) c[i-LA] = 0.0f;
  for( int I = std::max(1, LA-N+1); I <= std::min(M, LB); I++ ) {
    int MK = std::max(I, LA), K = MK-I+1;
    for( int J = MK; J <= std::min(N+I-1, LB); J++ ) { c[J-LA] = B[I-1]*A[K-1] + c[J-LA]; K++; }
  }
}
// MTRAPZ (DF = 1)
static void mtrapz( float* x, int n, float fl1, float fl2, float fh1, float fh2 ) {
  float rf1 = 1.0f/std::max( fl2-fl1, 1.0f ), rf2 = 1.0f/std::max( fh2-fh1, 1.0f );
  for( int i = 1; i <= n; i++ ) {
    float f = (float)(i-1);
    float t1 = std::min( std::max( rf1*(f-fl1), 0.0f ), 1.0f );
    float t2 = std::min( std::max( rf2*(fh2-f), 0.0f ), 1.0f );
    x[i-1] = x[i-1]*(t1+t2-1.0f);
  }
}
// LEVI: symmetric Toeplitz solve (A(1..N), Y destroyed, X result), N > 0
static int levi( float* A, int N, float* Y, float* X ) {
  int ier = 0;
  if( A[0] == 0.0f ) return 1;
  X[0] = Y[0]/A[0];
  Y[0] = A[0];
  for( int now = 2; now <= N; now++ ) {
    double temp = A[now-1], dot = temp, hold = Y[now-1];
    hold = hold - temp*X[0];
    for( int i = 2, k = now; i < now; i++ ) { k--; temp = A[k-1]; dot -= temp*Y[i-1]; hold -= temp*X[i-1]; }
    double bet = -dot/Y[0];
    if( fabs(bet) >= 1.0 ) ier = 1;
    Y[0] = (float)( Y[0] + bet*dot );
    Y[now-1] = (float)(-bet);
    for( int i = 2, k = now; ; ) { k--; if( i > k ) break; double d = Y[i-1], t = Y[k-1];
      Y[i-1] = (float)( d + bet*t ); Y[k-1] = (float)( t + bet*d ); i++; }
    bet = hold/Y[0];
    X[now-1] = (float)bet;
    for( int i = 2, k = now; i <= now; i++ ) { k--; X[k-1] = (float)( X[k-1] - bet*Y[i-1] ); }
  }
  return ier;
}

// design of one filter from z (window, length liw) and p; returns filter in f
static void design( VariableStruct* v, float const* z, float const* p, int liw, float* f ) {
  int nflt = v->nflt, nflth = -v->iprd, nfltm = nflt-1;
  std::vector<float> r( nflt ), y( nflt );
  correlg( z, liw, z, liw, &r[0], liw, liw+nfltm );
  correlg( p, liw, z, liw, &y[0], liw-nflth, liw-nflth+nfltm );
  if( v->pdiag != 0.0f ) r[0] = r[0]*(1.0f + v->pdiag/100.0f);
  for( int i = 0; i < nflt; i++ ) f[i] = 0.0f;
  if( levi( &r[0], nflt, &y[0], f ) ) v->ntrIer++;
  // store
  if( v->nfilt == 0 ) v->ffirst.assign( f, f+nflt );
  for( int i = 0; i < nflt; i++ ) v->fsum[i] += f[i];
  v->nfilt++;
}

// LSQSUB0 (whole trace) or LSQSUB (windows): fz = F*Z merged, length nt
static void match( VariableStruct* v, float const* zin, float const* pin, int nt, float* fz ) {
  int nflt = v->nflt, nflth = -v->iprd;
  std::vector<float> f( nflt );
  if( v->liw == 0 ) {                                           // LSQSUB0
    design( v, zin, pin, nt, &f[0] );
    std::vector<float> c( nt + nflt + 4, 0.0f );
    if( nflth >= 0 ) convg( &f[0], nflt, zin, nt, &c[0], nflth+1, nflth+nt );
    else             convg( &f[0], nflt, zin, nt, &c[-nflth], 1, nflth+nt );
    for( int k = 0; k < nt; k++ ) fz[k] = c[k];
    return;
  }
  // LSQSUB: 1-based work array s: z at s(1), p at s(1+NT2), conv at s(.+NT4), sum at s(.+NT6)
  int NT2 = 2*nt, NT4 = 4*nt, NT6 = 6*nt;
  std::vector<float> s( 8*nt + 4*nflt + 16, 0.0f );
  float* S = &s[0] - 1;                                          // S[1] = s[0]
  for( int k = 1; k <= nt; k++ ) { S[k] = zin[k-1]; S[k+NT2] = pin[k-1]; }
  int LIW = v->liw, LOW = v->low, LOV = v->lov + 1;
  int LIWH = LIW/2, LOWH = LOW/2, LIOWH = LIWH - LOWH;
  // number of windows
  int IST = 1, IEND = IST+LIW-1, JST = IST+LIOWH, JEND = JST+LOW-1, JJ = 0;
  for( int i = 1; i <= 500; i++ ) {
    JJ++;
    if( IEND == nt ) break;
    JST = JEND-LOV; JEND = JST+LOW-1; IST = JST-LIOWH; IEND = IST+LIW-1;
    if( JEND > nt-LOV ) break;
    if( IEND > nt ) IEND = nt;
  }
  IST = 1; IEND = IST+LIW-1; JST = IST+LIOWH; JEND = JST+LOW-1;
  int LEND = 0;
  for( int J = 1; J <= JJ; J++ ) {
    int liwUse = std::min( LIW, nt - IST + 1 );                 // SEASEIS: stay inside the trace
    design( v, &S[IST], &S[IST+NT2], liwUse, &f[0] );
    if( J == 1 ) {
      if( nflth >= 0 ) convg( &f[0], nflt, &S[IST], liwUse, &S[IST+NT4], nflth+1, LIOWH+nflth+LOW );
      else             convg( &f[0], nflt, &S[IST], liwUse, &S[IST+NT4-nflth], 1, LIOWH+nflth+LOW );
      float FH2 = JEND-IST, FH1 = JEND-IST-LOV;
      mtrapz( &S[IST+NT4], JEND, 0.f, 0.f, FH1, FH2 );
      for( int k = 0; k < JEND; k++ ) S[IST+NT6+k] += S[IST+NT4+k];
    }
    else if( J == JJ ) {
      if( LIOWH+nflth >= 0 ) convg( &f[0], nflt, &S[IST], liwUse, &S[JST+NT4], LIOWH+nflth+1, LIOWH+nflth+LOW+LEND );
      else convg( &f[0], nflt, &S[IST], liwUse, &S[JST+NT4-(LIOWH+nflth)], 1, LIOWH+nflth+LOW+LEND );
      float FL2 = LOV, FH2 = IEND-JST, FH1 = IEND-JST;
      mtrapz( &S[JST+NT4], LOW, 0.f, FL2, FH1, FH2 );
      for( int k = 0; k < LOW+LEND; k++ ) S[JST+NT6+k] += S[JST+NT4+k];
    }
    else {
      if( LIOWH+nflth >= 0 ) convg( &f[0], nflt, &S[IST], liwUse, &S[JST+NT4], LIOWH+nflth+1, LIOWH+nflth+LOW );
      else convg( &f[0], nflt, &S[IST], liwUse, &S[JST+NT4-(LIOWH+nflth)], 1, LIOWH+nflth+LOW );
      float FL2 = LOV, FH2 = JEND-JST, FH1 = JEND-JST-LOV;
      mtrapz( &S[JST+NT4], LOW, 0.f, FL2, FH1, FH2 );
      for( int k = 0; k < LOW; k++ ) S[JST+NT6+k] += S[JST+NT4+k];
    }
    JST = JEND-LOV; JEND = JST+LOW-1; IST = JST-LIOWH; IEND = IST+LIW-1;
    if( J == JJ-1 ) { IEND = nt; LIW = IEND-IST+1; LEND = IEND-JEND; }
  }
  for( int k = 1; k <= nt; k++ ) fz[k-1] = S[NT6+k];
}

//---------------------------------------------------------------- init
void init_mod_spk_lsqmatch_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* v = new VariableStruct();
  edef->setVariables( v );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  v->reader = NULL; v->hdef2 = NULL; v->hdrBlock = NULL;

  // defaults as elsqmatch.f
  float type = 1., len = 1., pdiag = 2., pred = 0., tinwl = 0., toutwl = 0., toutov = 0., t0 = 0., vel = 0., x;
  if( param->exists("type") )   param->getFloat("type", &type);
  if( param->exists("len") )    param->getFloat("len", &len);
  if( param->exists("pdiag") )  param->getFloat("pdiag", &pdiag);
  if( param->exists("pred") )   param->getFloat("pred", &pred);
  if( param->exists("tinwl") )  param->getFloat("tinwl", &tinwl);
  if( param->exists("toutwl") ) param->getFloat("toutwl", &toutwl);
  if( param->exists("toutov") ) param->getFloat("toutov", &toutov);
  if( param->exists("t0") )     param->getFloat("t0", &t0);
  if( param->exists("vel") )    param->getFloat("vel", &vel);
  const char* xpar[] = { "xinwl", "xoutwl", "xoutov" };
  for( int i = 0; i < 3; i++ ) if( param->exists(xpar[i]) ) { param->getFloat(xpar[i], &x); if( x != 0. )
    writer->error("'%s' (X windows, multi-trace) is not implemented yet", xpar[i]); }
  v->type = (int)type;
  if( v->type == 2 || v->type == 4 || v->type == 6 ) writer->error("TYPE %d (multi-trace) is not implemented yet: use TYPE 1, 3 or 5", v->type);
  if( v->type != 1 && v->type != 3 && v->type != 5 ) writer->error("TYPE must be 1, 3 or 5");

  v->nt = shdr->numSamples; v->sr = shdr->sampleInt;
  v->nflt = (int)len;
  if( v->nflt >= v->nt/2 ) v->nflt = v->nt/2 - 1;
  if( v->nflt < 1 ) writer->error("LEN must be >= 1");
  v->iprd = (int)pred;
  v->pdiag = pdiag;
  v->liw = (int)lround( tinwl/v->sr );
  v->low = (int)lround( toutwl/v->sr );
  v->lov = (int)lround( toutov/v->sr );
  if( v->liw > 0 && v->low <= 0 ) writer->error("TOUTWL must be > 0 when TINWL is given");
  if( v->liw > 0 && v->low > v->liw ) writer->error("TOUTWL must not be larger than TINWL");
  if( v->liw > v->nt ) v->liw = v->nt;
  if( v->low > v->nt ) v->low = v->nt;
  v->t0 = t0; v->vel = vel;
  v->idOff = -1;
  if( t0 > 0. && vel > 0. ) {
    string name = hdef->headerExists("spk_rang") ? "spk_rang" : "offset";
    if( param->exists("hdr_offset") ) param->getString("hdr_offset", &name);
    if( !hdef->headerExists(name) ) writer->error("Offset header '%s' not found (needed for T0/VEL): set hdr_offset", name.c_str());
    v->idOff = hdef->headerIndex( name );
  }

  // ONLINE2
  string fname;
  if( !param->exists("filename") ) writer->error("Parameter 'filename' (ONLINE2, data to be matched) is required");
  param->getString( "filename", &fname );
  v->hdef2 = new csTraceHeaderDef( hdef );
  v->reader = new csSeismicReader( fname, false, 1 );
  csSuperHeader shdr2; int blk = 0;
  if( !v->reader->readFileHeader( &shdr2, v->hdef2, &blk ) ) writer->error("Cannot read file header of '%s'", fname.c_str());
  v->hdef2->resetByteLocation();
  v->hdrBlock = new char[ blk + 8 ];
  v->nt2 = shdr2.numSamples;
  if( fabs( shdr2.sampleInt - v->sr ) > 1e-4 ) writer->error("ONLINE2 '%s' has sample interval %g, ONLINE1 %g", fname.c_str(), shdr2.sampleInt, v->sr);
  if( v->nt2 != v->nt ) writer->warning("ONLINE2 has %d samples, ONLINE1 %d: missing samples are zero", v->nt2, v->nt);
  v->z.assign( std::max( v->nt, v->nt2 ) + 8, 0.0f );

  v->filterFile = "FILTER";
  if( param->exists("filter_file") ) param->getString("filter_file", &v->filterFile);
  string fm = "first";
  if( param->exists("filter") ) param->getString("filter", &fm);
  if( fm != "first" && fm != "mean" ) writer->error("filter must be 'first' or 'mean'");
  v->filterMean = ( fm == "mean" );
  v->fsum.assign( v->nflt, 0.0 ); v->nfilt = 0; v->ntrZeroZ = 0; v->ntrIer = 0;

  writer->line("  ONLINE2 (matched to the flow): %s  (%d traces, %d samples)", fname.c_str(), v->reader->numTraces(), v->nt2);
  writer->line("  TYPE %d: output = %s", v->type, v->type == 1 ? "ONLINE1 - F*ONLINE2" : v->type == 3 ? "ONLINE1 + F*ONLINE2" : "F*ONLINE2");
  writer->line("  LEN %d  PRED %d  PDIAG %g %%", v->nflt, v->iprd, v->pdiag);
  if( v->liw == 0 ) writer->line("  Whole trace design (TINWL = 0)");
  else writer->line("  Windows [samples]: in %d  out %d  overlap %d", v->liw, v->low, v->lov);
  if( v->idOff >= 0 ) writer->line("  Moveout shift: T0 %g ms VEL %g m/s", v->t0, v->vel);
  writer->line("  FILTER file: %s (%s filter)", v->filterFile.c_str(), v->filterMean ? "mean" : "first");
}

//---------------------------------------------------------------- exec
void exec_mod_spk_lsqmatch_( csTraceGather* traceGather, int* port, int* numTrcToKeep,
                             csExecPhaseEnv* env, csLogWriter* writer )
{
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csTrace* trace = traceGather->trace(0);
  float* p = trace->getTraceSamples();
  int nt = v->nt;
  for( size_t k = 0; k < v->z.size(); k++ ) v->z[k] = 0.0f;
  if( !v->reader->readTrace( &v->z[0], v->hdrBlock, v->nt2 ) ) {
    writer->error("ONLINE2 has fewer traces than the flow (ONLINE1)");
  }
  float* z = &v->z[0];
  float dmax = 0.0f;
  for( int j = 0; j < nt; j++ ) dmax = std::max( dmax, fabsf( p[j] ) );
  if( dmax < .001f ) { for( int j = 0; j < nt; j++ ) z[j] = 0.0f; v->ntrZeroZ++; }

  int ishft = 0;
  if( v->idOff >= 0 ) {
    double rng = trace->getTraceHeader()->doubleValue( v->idOff );
    double tt = sqrt( v->t0*v->t0 + pow( 1000.*rng/v->vel, 2 ) );
    ishft = (int)lround( (tt - v->t0)/v->sr );
    if( ishft < 0 ) ishft = 0;
    if( ishft > nt ) ishft = nt;
  }
  std::vector<float> ps( nt, 0.0f ), zs( nt, 0.0f ), fz( nt, 0.0f );
  for( int j = 0; j < nt-ishft; j++ ) { ps[j] = p[j+ishft]; zs[j] = z[j+ishft]; }
  match( v, &zs[0], &ps[0], nt, &fz[0] );

  std::vector<float> outs( nt );
  for( int j = 0; j < nt; j++ ) outs[j] = ( v->type == 3 ) ? ps[j] + fz[j] : ps[j] - fz[j];
  std::vector<float> out( nt, 0.0f );
  for( int j = 0; j < nt-ishft; j++ ) out[j+ishft] = outs[j];
  if( v->type == 5 ) for( int j = 0; j < nt; j++ ) out[j] = p[j] - out[j];    // ONLINE1 - (ONLINE1 - F*ONLINE2)
  memcpy( p, &out[0], nt*sizeof(float) );
}

//---------------------------------------------------------------- params
void params_mod_spk_lsqmatch_( csParamDef* pdef ) {
  pdef->setModule( "SPK_LSQMATCH", "Least-squares match filter (SEISPAK lsqmatch)",
                   "Match the traces of 'filename' (ONLINE2) to the flow (ONLINE1) trace by trace with a Wiener filter "
                   "(Levinson), on the whole trace or in sliding time windows. TYPE 1: ONLINE1 - F*ONLINE2, TYPE 3: ONLINE1 + F*ONLINE2, "
                   "TYPE 5: F*ONLINE2. The filter is written to FILTER (for SPK_MFILT2)." );
  pdef->addParam( "filename", "ONLINE2: file with the data to be matched (.cseis), same trace order as the flow", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addParam( "type", "1: ONLINE1-F*ONLINE2, 3: ONLINE1+F*ONLINE2, 5: F*ONLINE2 (trace-wise)", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "TYPE" );
  pdef->addParam( "len", "Filter length [samples]", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "LEN" );
  pdef->addParam( "pred", "Prediction distance [samples] (-(LEN/2) for a match filter)", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "PRED" );
  pdef->addParam( "pdiag", "Percent added to the zero lag of the autocorrelation", NUM_VALUES_FIXED );
  pdef->addValue( "2", VALTYPE_NUMBER, "PDIAG [%]" );
  pdef->addParam( "tinwl", "Input (design) window length [ms] (0: whole trace)", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "TINWL [ms]" );
  pdef->addParam( "toutwl", "Output window length [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "TOUTWL [ms]" );
  pdef->addParam( "toutov", "Output window overlap [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "TOUTOV [ms]" );
  pdef->addParam( "t0", "T0 for the moveout shift [ms] (with VEL)", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "T0 [ms]" );
  pdef->addParam( "vel", "Velocity for the moveout shift [m/s]", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "VEL [m/s]" );
  pdef->addParam( "hdr_offset", "Offset header (RANG) for the moveout shift", NUM_VALUES_FIXED, "Default: spk_rang, else offset" );
  pdef->addValue( "spk_rang", VALTYPE_STRING, "Header name" );
  pdef->addParam( "filter_file", "Output filter file (Fortran unformatted, as SEISPAK FILTER)", NUM_VALUES_FIXED );
  pdef->addValue( "FILTER", VALTYPE_STRING, "File name" );
  pdef->addParam( "filter", "Filter written to filter_file", NUM_VALUES_FIXED );
  pdef->addValue( "first", VALTYPE_OPTION );
  pdef->addOption( "first", "First filter computed (as xlsqmatch)" );
  pdef->addOption( "mean", "Average of all filters (traces and windows)" );
  pdef->addParam( "xinwl", "X windows (multi-trace): not implemented", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "" );
  pdef->addParam( "xoutwl", "X windows (multi-trace): not implemented", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "" );
  pdef->addParam( "xoutov", "X windows (multi-trace): not implemented", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "" );
}

bool start_exec_mod_spk_lsqmatch_( csExecPhaseEnv* env, csLogWriter* writer ) { return true; }

static void write_rec( FILE* f, void const* data, int nbytes ) {
  int n = nbytes; fwrite( &n, 4, 1, f ); fwrite( data, 1, nbytes, f ); fwrite( &n, 4, 1, f );
}

void cleanup_mod_spk_lsqmatch_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( v == NULL ) return;
  if( v->nfilt > 0 ) {
    std::vector<float> f( v->nflt );
    for( int i = 0; i < v->nflt; i++ ) f[i] = v->filterMean ? (float)( v->fsum[i]/v->nfilt ) : v->ffirst[i];
    FILE* fo = fopen( v->filterFile.c_str(), "wb" );
    if( fo ) {
      struct { int nf, ntr, nrec; float sr; } h = { v->nflt, v->nfilt, 1, 1.0f };
      write_rec( fo, &h, 16 );
      write_rec( fo, &f[0], 4*v->nflt );
      fclose( fo );
      writer->line("  FILTER '%s' written: %d coefficients (%s of %d filters)", v->filterFile.c_str(), v->nflt,
                   v->filterMean ? "mean" : "first", v->nfilt);
    }
    else writer->warning("Cannot write filter file '%s'", v->filterFile.c_str());
    if( v->ntrZeroZ > 0 ) writer->line("  %d zero ONLINE1 traces (ONLINE2 zeroed)", v->ntrZeroZ);
    if( v->ntrIer > 0 ) writer->line("  Levinson: %d filters with |reflection coefficient| >= 1 or zero autocorrelation", v->ntrIer);
  }
  delete v->reader; delete v->hdef2; delete [] v->hdrBlock;
  delete v; v = NULL;
}

extern "C" void _params_mod_spk_lsqmatch_( csParamDef* pdef ) { params_mod_spk_lsqmatch_( pdef ); }
extern "C" void _init_mod_spk_lsqmatch_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spk_lsqmatch_( param, env, writer );
}
extern "C" bool _start_exec_mod_spk_lsqmatch_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spk_lsqmatch_( env, writer );
}
extern "C" void _exec_mod_spk_lsqmatch_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_lsqmatch_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_lsqmatch_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spk_lsqmatch_( env, writer );
}
