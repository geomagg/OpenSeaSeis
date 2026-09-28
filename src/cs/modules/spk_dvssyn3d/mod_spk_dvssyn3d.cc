/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module SPK_DVSSYN3D: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

#include "csStandardHeaders.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_DVSSYN3D
 *
 * Input module: 2D/3D zero-offset synthetic (exploding reflector) from SEISPAK
 * program DVSSYN. Reflectivity and velocities are built from the velocity model
 * (VXT horizon list or SEISPAK ONL v(z) file) by XMODELI2, and the wavefield is
 * upward-continued with the DVSSYN operators (phase shift, 45-degree x-w, PSPI,
 * multi-velocity PSPI; 3D: dwn3d, phs2d, mvphs2d). Parameters follow edvssyn/pdvssyn.
 */

extern "C" {
  void spksyn_( float* vel, int* nvelw, float* afit, int* nafit, int* ivt, int* nx, int* ny, int* nt, float* sr,
                int* nnz, int* itzr, float* dx, float* dy, float* ddz, float* apx, float* alph1, float* oty,
                float* frq1, float* frq2, int* irfc, float* out, int* info, int* ierr );
  void spksynvxt_( float* vel, float* dx, int* nxv, int* np, float* afit, float* ee, float* x, float* bfit, float* xx );
}

namespace mod_spk_dvssyn3d {
  struct VariableStruct {
    int nx, ny, nt;
    std::vector<float> out;
    int traceCounter;
    bool atEOF;
    bool computed;
    // parameters
    float sr, dx, dy, ddz, apx, alph1, oty, frq1, frq2;
    int nnz, itzr, irfc, ivt;
    std::vector<float> vel;
    std::vector<float> afit;
    int hdrId_trcno, hdrId_row, hdrId_col, hdrId_cmp;
  };

  static int const BASE = 12;
  static inline float getF( char const* p, bool sw ) {
    char b[4]; memcpy(b,p,4); if(sw){ char c=b[0];b[0]=b[3];b[3]=c; c=b[1];b[1]=b[2];b[2]=c; }
    float v; memcpy(&v,b,4); return v;
  }
  static inline int getI( char const* p, bool sw ) {
    char b[4]; memcpy(b,p,4); if(sw){ char c=b[0];b[0]=b[3];b[3]=c; c=b[1];b[1]=b[2];b[2]=c; }
    int v; memcpy(&v,b,4); return v;
  }
  static inline short getS( char const* p, bool sw ) {
    char b[2]; memcpy(b,p,2); if(sw){ char c=b[0];b[0]=b[1];b[1]=c; }
    short v; memcpy(&v,b,2); return v;
  }
  // Reads a SEISPAK file: all traces, nt samples, sample interval, HMT (traces/record), HMR
  std::string readSeispak( std::string const& fname, int* hmtOut, int* hmrOut, int* ntOut, float* smpOut,
                           std::vector<float>* data ) {
    FILE* fp = fopen( fname.c_str(), "rb" );
    if( fp == NULL ) return "cannot open file";
    std::vector<char> raw;
    char buf[65536];
    size_t n;
    while( (n = fread(buf,1,sizeof(buf),fp)) > 0 ) raw.insert( raw.end(), buf, buf+n );
    fclose(fp);
    if( raw.size() < 200 || strncmp(&raw[0],"ICE",3) != 0 ) return "not a SEISPAK file (no 'ICE' signature)";
    char const* d = &raw[BASE];
    bool sw = false;
    int lh = getI(d+4,false);
    if( !(lh >= 0 && lh < 10000) ) { sw = true; lh = getI(d+4,true); }
    int ifw0 = getI(d,sw), nt = getI(d+8,sw), nbyte = getI(d+32,sw);
    float smp = getF(d+12,sw), hmt = getF(d+16,sw), hmr = getF(d+20,sw);
    if( nbyte == 0 ) nbyte = 4;
    int ihmt = (int)(hmt+0.5f), ihmr = (int)(hmr+0.5f);
    int ntr = ihmt*ihmr;
    if( nt <= 0 || ntr <= 0 || ifw0 < 21 ) return "implausible header directory";
    long long rec;
    if( nbyte == 4 ) rec = 4LL*(lh+nt);
    else if( nbyte == 2 ) rec = 4LL*(lh + nt/2 + ((nt%2==0)?1:2));
    else return "only NBYTE 2 and 4 supported for velocity files";
    long long first = BASE + 4LL*(ifw0-1);
    if( first + rec*ntr > (long long)raw.size() ) return "file truncated";
    data->assign( (size_t)ntr*nt, 0.0f );
    for( int k = 0; k < ntr; k++ ) {
      char const* p = &raw[first + rec*k] + 4*lh;
      float* out = &(*data)[(size_t)k*nt];
      if( nbyte == 4 ) {
        for( int i = 0; i < nt; i++ ) out[i] = getF(p+4*i,sw);
      }
      else {
        float fac = getF(p,sw);
        char const* s = p + 4 + ((nt%2) ? 2 : 0);
        for( int i = 0; i < nt; i++ ) out[i] = (fac > 1e-30f) ? getS(s+2*i,sw)/fac : 0.0f;
      }
    }
    *hmtOut = ihmt; *hmrOut = ihmr; *ntOut = nt; *smpOut = smp;
    return "";
  }
}
using namespace mod_spk_dvssyn3d;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_spk_dvssyn3d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  edef->setExecType( EXEC_TYPE_INPUT );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->traceCounter = 0;
  vars->atEOF = false;
  vars->computed = false;
  // Presets as in edvssyn.f
  vars->nx = 0; vars->ny = 1; vars->nt = 0; vars->nnz = 0;
  vars->itzr = 5; vars->sr = 4.0f; vars->dx = 0.0f; vars->dy = 0.0f; vars->ddz = 0.0f;
  vars->apx = 4.0f; vars->alph1 = 1.0f; vars->oty = 0.0f; vars->frq1 = 10.0f; vars->frq2 = 60.0f;
  vars->irfc = 1; vars->ivt = 0;

  std::string text;
  if( param->exists("nx") )    param->getInt( "nx", &vars->nx );
  if( param->exists("ny") )    param->getInt( "ny", &vars->ny );
  if( param->exists("nt") )    param->getInt( "nt", &vars->nt );
  if( param->exists("nnz") )   param->getInt( "nnz", &vars->nnz );
  if( param->exists("itzr") )  param->getInt( "itzr", &vars->itzr );
  if( param->exists("sr") )    param->getFloat( "sr", &vars->sr );
  if( param->exists("dx") )    param->getFloat( "dx", &vars->dx );
  if( param->exists("dy") )    param->getFloat( "dy", &vars->dy );
  if( param->exists("ddz") )   param->getFloat( "ddz", &vars->ddz );
  if( param->exists("apx") )   param->getFloat( "apx", &vars->apx );
  if( param->exists("alph1") ) param->getFloat( "alph1", &vars->alph1 );
  if( param->exists("oty") )   param->getFloat( "oty", &vars->oty );
  if( param->exists("frq1") )  param->getFloat( "frq1", &vars->frq1 );
  if( param->exists("frq2") )  param->getFloat( "frq2", &vars->frq2 );
  if( param->exists("rfc") ) {
    param->getString( "rfc", &text );
    if( !text.compare("yes") || !text.compare("YES") ) vars->irfc = 1;
    else if( !text.compare("no") || !text.compare("NO") ) vars->irfc = 0;
    else writer->error("Unknown option for 'rfc': %s", text.c_str());
  }
  if( vars->ny < 1 ) vars->ny = 1;

  // Required parameters (pdvssyn)
  if( vars->nx <= 0 )  writer->error("Required parameter NX (depth points in x) not given");
  if( vars->nnz <= 0 ) writer->error("Required parameter NNZ (depth samples) not given");
  if( vars->nt <= 0 )  writer->error("Required parameter NT (time samples of output) not given");
  if( vars->dx <= 0 )  writer->error("Required parameter DX (depth point spacing) not given");
  if( vars->ddz <= 0 ) writer->error("Required parameter DDZ (depth sample rate) not given");
  if( vars->sr <= 0 )  writer->error("SR must be > 0");
  if( vars->itzr < 1 || vars->nnz < vars->itzr ) writer->error("ITZR must be >= 1 and <= NNZ");
  if( vars->frq2 <= vars->frq1 ) writer->error("FRQ2 must be larger than FRQ1");
  if( vars->alph1 < 1.0f ) writer->error("ALPH1 must be >= 1");
  if( (int)lrintf(vars->alph1) * vars->nt > 9000 ) writer->error("NT*ALPH1 must be <= 9000 (array limit of XMODELI2)");
  float a = vars->apx;
  bool is2D = ( vars->ny == 1 );
  bool apxOK = is2D ? ( a == 1.0f || a == 4.0f || a == 5.0f || a == 6.0f || a == 7.0f || a == 10.0f || a >= 11.0f )
                    : ( a == 4.0f || a == 5.0f || a == 6.0f || a == 7.0f || a >= 11.0f );
  if( !apxOK ) writer->error("APX = %g not supported for %s modelling", a, is2D ? "2D" : "3D");
  if( a >= 11.0f && (int)lrintf(a-10.0f) < 2 ) writer->error("APX >= 11 requires at least 2 velocities (APX >= 12)");
  if( !is2D && vars->dy <= 0.0f ) writer->warning("DY not given: DY = DX used");

  //---------------------------------------------
  // Velocity model
  bool hasVxt = param->exists("vxt_file"), hasVel = param->exists("vel_file");
  if( hasVxt == hasVel ) writer->error("Specify one velocity source: 'vxt_file' (VXT list) or 'vel_file' (SEISPAK ONL v(z) file)");
  if( hasVxt ) {
    std::string fname;
    param->getString( "vxt_file", &fname );
    FILE* fp = fopen( fname.c_str(), "r" );
    if( fp == NULL ) writer->error("Cannot open VXT file '%s'", fname.c_str());
    std::vector<std::string> tok;
    char word[256];
    while( fscanf( fp, "%255s", word ) == 1 ) tok.push_back( word );
    fclose( fp );
    size_t start = 0;
    for( size_t i = 0; i < tok.size(); i++ ) {
      std::string t = tok[i];
      for( size_t c = 0; c < t.size(); c++ ) t[c] = (char)tolower(t[c]);
      if( t == "vxt" || t == "vel" ) { start = i+1; break; }
    }
    std::vector<float> v;
    for( size_t i = start; i < tok.size(); i++ ) {
      char* end = NULL;
      float f = strtof( tok[i].c_str(), &end );
      if( end == tok[i].c_str() || *end != 0 ) break;
      v.push_back( f );
    }
    if( v.size() < 5 || v[0] <= 0.0f ) writer->error("VXT file '%s': no velocity list found (first value must be the VXT code > 0)", fname.c_str());
    // Split into y-line groups: [code] [negative limits] y1 <horizons...> 0 y2 <horizons...> 0 ...
    size_t pos = 1;
    std::vector<float> lead;
    while( pos < v.size() && v[pos] < 0.0f ) lead.push_back( v[pos++] );
    std::vector<float> ys;
    std::vector< std::vector<float> > bodies;
    while( pos < v.size() ) {
      float y = v[pos++];
      size_t s0 = pos;
      while( pos+1 < v.size() && !( v[pos] <= 0.5f && v[pos+1] <= 0.5f ) ) pos++;
      bodies.push_back( std::vector<float>( v.begin()+s0, v.begin()+std::min(pos+1, v.size()) ) );
      ys.push_back( y );
      pos += 2;
    }
    int ngroups = (int)ys.size();
    if( ngroups == 0 ) writer->error("VXT file '%s': no y-line found", fname.c_str());
    for( int g = 1; g < ngroups; g++ ) if( ys[g] <= ys[g-1] ) writer->error("VXT file: y-lines must increase (y=%g after y=%g)", ys[g], ys[g-1]);
    // Dense functions along x for each y-line (original XVEVENTS)
    std::vector< std::vector<float> > gafit( ngroups );
    std::vector<int> gnx( ngroups ), gnp( ngroups );
    for( int g = 0; g < ngroups; g++ ) {
      std::vector<float> list;
      list.push_back( v[0] );
      list.insert( list.end(), lead.begin(), lead.end() );
      list.push_back( 1.0f );
      list.insert( list.end(), bodies[g].begin(), bodies[g].end() );
      list.push_back( 0.0f );
      size_t m = list.size() + 60000;
      list.resize( m, 0.0f );
      std::vector<float> ee( m ), x( 20*m + 100000 ), xx( 10001 );
      size_t na = 3 + (size_t)10001*203;
      gafit[g].assign( na, 0.0f );
      std::vector<float> bfit( na );
      float ddd = vars->dx;
      spksynvxt_( &list[0], &ddd, &gnx[g], &gnp[g], &gafit[g][0], &ee[0], &x[0], &bfit[0], &xx[0] );
      if( gnx[g] <= 0 || gnp[g] <= 0 || gnx[g] > 10000 ) writer->error("VXT file: could not build velocity functions for y-line %g", ys[g]);
      if( gnx[g] != gnx[0] ) writer->error("VXT file: y-line %g has %d DPN, y-line %g has %d. All y-lines must span the same DPN range.", ys[g], gnx[g], ys[0], gnx[0]);
    }
    int nxv = gnx[0], np = 0;
    for( int g = 0; g < ngroups; g++ ) np = std::max( np, gnp[g] );
    int nyv = std::max( 1, vars->ny );
    int nw = 2*np+3;
    vars->afit.assign( 3 + (size_t)nxv*nyv*nw + 1000, 0.0f );
    vars->afit[0] = (float)nxv; vars->afit[1] = (float)nyv; vars->afit[2] = (float)np;
    for( int iy = 1; iy <= nyv; iy++ ) {
      // bracketing y-lines (model y index = y value of the VXT line)
      int g0 = 0, g1 = 0;
      float w = 0.0f;
      if( ngroups > 1 ) {
        if( iy <= ys[0] ) { g0 = g1 = 0; }
        else if( iy >= ys[ngroups-1] ) { g0 = g1 = ngroups-1; }
        else {
          g1 = 1;
          while( g1 < ngroups-1 && ys[g1] < iy ) g1++;
          g0 = g1-1;
          w = ( iy - ys[g0] ) / ( ys[g1] - ys[g0] );
        }
      }
      if( gnp[g0] != gnp[g1] && w > 0.0f ) { if( w >= 0.5f ) g0 = g1; w = 0.0f; }
      int nw0 = 2*gnp[g0]+3, nw1 = 2*gnp[g1]+3;
      for( int ix = 1; ix <= nxv; ix++ ) {
        float const* r0 = &gafit[g0][3 + (size_t)(ix-1)*nw0];
        float const* r1 = &gafit[g1][3 + (size_t)(ix-1)*nw1];
        float* r = &vars->afit[3 + ((size_t)(iy-1)*nxv + ix-1)*nw];
        r[0] = (float)iy; r[1] = r0[1];
        for( int k = 0; k < 2*gnp[g0]; k++ ) r[2+k] = (1.0f-w)*r0[2+k] + w*r1[2+k];
      }
    }
    vars->vel.assign( 60000, 0.0f );
    vars->vel[0] = -3.0f;
    vars->ivt = 3;
    float ymax = ys.back(), dpnmax = (float)nxv; int nmMax = np;
    if( ngroups > 1 && vars->ny < ymax ) writer->warning("NY (%d) smaller than the last y-line of the VXT (%g)", vars->ny, ymax);
    if( ngroups > 1 && vars->ny == 1 ) writer->warning("VXT has %d y-lines but NY = 1: only y-line %g is used", ngroups, ys[0]);
    writer->line("  VXT file: %s  (code %g, %d y-lines, up to %d horizons, DPN up to %g)", fname.c_str(), v[0], ngroups, nmMax, dpnmax);
  }
  else {
    std::string fname;
    param->getString( "vel_file", &fname );
    int hmt, hmr, ntv; float dzv;
    std::vector<float> vdat;
    std::string err = readSeispak( fname, &hmt, &hmr, &ntv, &dzv, &vdat );
    if( !err.empty() ) writer->error("Velocity file '%s': %s", fname.c_str(), err.c_str());
    // As ONLVELS: record = y line, trace = DPN. For a 2D file with one trace per record, records are the DPNs.
    int nxv = hmt, nyv = hmr;
    if( hmt == 1 && hmr > 1 && vars->ny == 1 ) { nxv = hmr; nyv = 1; }
    int np = ntv, nw = 2*np+3;
    vars->afit.assign( 3 + (size_t)nxv*nyv*nw + 1000, 0.0f );
    vars->afit[0] = (float)nxv; vars->afit[1] = (float)nyv; vars->afit[2] = (float)np;
    for( int iy = 1; iy <= nyv; iy++ ) {
      for( int ix = 1; ix <= nxv; ix++ ) {
        float* xx = &vars->afit[3 + ((size_t)(iy-1)*nxv + ix-1)*nw];
        float const* tr = &vdat[((size_t)(iy-1)*nxv + ix-1)*ntv];
        xx[0] = (float)iy; xx[1] = (float)ix;
        for( int j = 0; j < np; j++ ) { xx[2+2*j] = tr[j]; xx[3+2*j] = j*dzv; }
      }
    }
    vars->vel.assign( 60000, 0.0f );
    vars->vel[0] = -3.0f;
    vars->ivt = 3;
    writer->line("  Velocity file (ONL): %s  (%d DPN x %d y-lines, %d samples, depth step %f)", fname.c_str(), nxv, nyv, np, dzv);
  }

  //---------------------------------------------
  shdr->numSamples = vars->nt;
  shdr->sampleInt  = vars->sr;
  shdr->domain     = DOMAIN_XT;
  vars->hdrId_trcno = hdef->addStandardHeader( HDR_TRCNO.name );
  vars->hdrId_cmp   = hdef->addStandardHeader( HDR_CMP.name );
  vars->hdrId_row   = hdef->addStandardHeader( HDR_ROW.name );
  vars->hdrId_col   = hdef->addStandardHeader( HDR_COL.name );

  char const* apxText = ( a == 1.0f ) ? "phase shift" : ( a < 6.0f ) ? "45-degree x-w" : ( a < 10.0f ) ? "PSPI" :
                        ( a == 10.0f ) ? "table-driven convolution" : "multi-velocity PSPI";
  writer->line("  %s zero-offset synthetic: NX %d  NY %d  NT %d  SR %f ms", is2D ? "2D" : "3D", vars->nx, vars->ny, vars->nt, vars->sr);
  writer->line("  APX %g (%s)  ITZR %d  NNZ %d  DDZ %f  DX %f  DY %f  ALPH1 %g  OTY %g  FRQ1/FRQ2 %g %g  RFC %s",
               a, apxText, vars->itzr, vars->nnz, vars->ddz, vars->dx, vars->dy, vars->alph1, vars->oty, vars->frq1, vars->frq2,
               vars->irfc ? "YES" : "NO");
  writer->line("  Output: traces ordered by y-line (header 'row') and DPN (headers 'col' and 'cmp')");
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_spk_dvssyn3d_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( vars->atEOF ) {
    traceGather->freeAllTraces();
    return;
  }
  int ntr = vars->nx * vars->ny;
  if( !vars->computed ) {
    vars->out.assign( (size_t)ntr * vars->nt, 0.0f );
    int info[8] = {0,0,0,0,0,0,0,0};
    int ierr = 0;
    int nvelw = (int)vars->vel.size();
    int nafit = (int)vars->afit.size();
    writer->line("  Computing synthetic...");
    spksyn_( &vars->vel[0], &nvelw, &vars->afit[0], &nafit, &vars->ivt, &vars->nx, &vars->ny, &vars->nt, &vars->sr,
             &vars->nnz, &vars->itzr, &vars->dx, &vars->dy, &vars->ddz, &vars->apx, &vars->alph1, &vars->oty,
             &vars->frq1, &vars->frq2, &vars->irfc, &vars->out[0], info, &ierr );
    writer->line("  nx/ny (FFT) %d %d  nt1 %d  frequency samples %d-%d  depth steps %d  memory ~%d MB",
                 info[0], info[1], info[2], info[3], info[4], info[5], info[6]);
    if( ierr == 1 ) writer->error("No frequencies between FRQ1 and FRQ2");
    if( ierr == 2 ) writer->error("Not enough memory (%d MB)", info[6]);
    if( ierr == 3 ) writer->error("NT*ALPH1 too large");
    vars->afit.clear();
    vars->vel.clear();
    vars->computed = true;
  }
  int k = vars->traceCounter;
  if( k == ntr-1 ) vars->atEOF = true;
  csTrace* trace = traceGather->trace(0);
  memcpy( trace->getTraceSamples(), &vars->out[(size_t)k*vars->nt], vars->nt*sizeof(float) );
  csTraceHeader* trcHdr = trace->getTraceHeader();
  int iy = k / vars->nx + 1;
  int ix = k % vars->nx + 1;
  trcHdr->setIntValue( vars->hdrId_trcno, k+1 );
  trcHdr->setIntValue( vars->hdrId_row, iy );
  trcHdr->setIntValue( vars->hdrId_col, ix );
  trcHdr->setIntValue( vars->hdrId_cmp, ix );
  vars->traceCounter += 1;
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_spk_dvssyn3d_( csParamDef* pdef ) {
  pdef->setModule( "SPK_DVSSYN3D", "2D/3D zero-offset synthetic, exploding reflector (SEISPAK DVSSYN)",
                   "Input module. Builds reflectivity (velocity contrasts) and velocities from the model and upward-continues "
                   "the exploding-reflector wavefield with the DVSSYN operators. NY = 1 gives a 2D line. "
                   "Parameters and presets follow edvssyn/pdvssyn. Computation is done by the original DVSSYN Fortran routines." );
  pdef->addParam( "nx", "Depth points in x direction", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NX" );
  pdef->addParam( "ny", "Depth points in y direction (1 = 2D)", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "NY" );
  pdef->addParam( "nnz", "Depth samples of the depth model", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NNZ" );
  pdef->addParam( "nt", "Time samples of the output synthetic", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NT" );
  pdef->addParam( "sr", "Time sample rate of the output [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "4", VALTYPE_NUMBER, "SR" );
  pdef->addParam( "dx", "Depth point spacing in x (DDD)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "DX" );
  pdef->addParam( "dy", "Depth point spacing in y (3D)", NUM_VALUES_FIXED, "0 = DX" );
  pdef->addValue( "0", VALTYPE_NUMBER, "DY" );
  pdef->addParam( "ddz", "Depth sample rate", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "DDZ" );
  pdef->addParam( "itzr", "Depth samples per continuation step", NUM_VALUES_FIXED );
  pdef->addValue( "5", VALTYPE_NUMBER, "ITZR" );
  pdef->addParam( "apx", "Approximation", NUM_VALUES_FIXED );
  pdef->addValue( "4", VALTYPE_NUMBER, "1: phase shift (2D), 4/5: 45-degree x-w, 6/7: PSPI (3D: phs2d), 10: table-driven (2D), "
                  ">=11: multi-velocity PSPI with APX-10 velocities, >=13 adds absorbing boundaries" );
  pdef->addParam( "alph1", "Refinement factor for the depth model", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "ALPH1" );
  pdef->addParam( "oty", "Shift from x=0 in the model to DPN 1 (in DPNs)", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "OTY" );
  pdef->addParam( "frq1", "Lower frequency limit", NUM_VALUES_FIXED );
  pdef->addValue( "10", VALTYPE_NUMBER, "Hz" );
  pdef->addParam( "frq2", "Upper frequency limit", NUM_VALUES_FIXED );
  pdef->addValue( "60", VALTYPE_NUMBER, "Hz" );
  pdef->addParam( "rfc", "RFC flag (used by APX 10 only)", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Depth model" );
  pdef->addOption( "no", "Time model" );
  pdef->addParam( "vxt_file", "File with SEISPAK VXT velocity list (horizons), as in the jobdeck", NUM_VALUES_FIXED,
                  "Numbers after the keyword VXT (or the whole file): code, then for each y-line: y, horizons as triplets "
                  "(velocity, DPN, depth or time) each ended by a value <= 0.5, and 0 at the end of the line. "
                  "Each y-line is interpolated along x by the original XVEVENTS; between y-lines, depth and velocity of "
                  "each horizon are interpolated linearly in y." );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addParam( "vel_file", "SEISPAK ONL file with interval velocity v(z)", NUM_VALUES_FIXED,
                  "One trace per DPN, one record per y-line; sample interval = depth step" );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
}

//************************************************************************************************
bool start_exec_mod_spk_dvssyn3d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

void cleanup_mod_spk_dvssyn3d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  delete vars; vars = NULL;
}

extern "C" void _params_mod_spk_dvssyn3d_( csParamDef* pdef ) {
  params_mod_spk_dvssyn3d_( pdef );
}
extern "C" void _init_mod_spk_dvssyn3d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spk_dvssyn3d_( param, env, writer );
}
extern "C" bool _start_exec_mod_spk_dvssyn3d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spk_dvssyn3d_( env, writer );
}
extern "C" void _exec_mod_spk_dvssyn3d_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_dvssyn3d_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_dvssyn3d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spk_dvssyn3d_( env, writer );
}
