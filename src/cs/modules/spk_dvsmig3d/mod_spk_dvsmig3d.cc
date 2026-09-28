/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module SPK_DVSMIG3D: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/types.h>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_DVSMIG3D
 *
 * 3D post-stack wave-equation migration from SEISPAK program DVSMIG (path ny > 1):
 * split 45-degree x-y-w finite difference (APX 4/5), phase shift plus split-step
 * correction (APX 6/7) and multi-velocity PSPI (APX >= 11), in time (RFC NO)
 * or depth (RFC YES). The whole volume is collected, organised in y-lines,
 * migrated by the Fortran code in fortran/, and output with the input headers.
 */

extern "C" {
  void spkmig3d_( float* d, int* nxsave, int* nysave, int* nt, float* sr, float* vxz, int* nz, int* itzr, int* nnz,
                  float* ddz, float* dx, float* dy, float* apx, int* irfc, float* frq1, float* frq2,
                  float* out, int* info, int* ierr );
  void spkvcol_( float* xx, int* np, int* nz, int* itzr, float* ddz, float* dt, int* irfc, float* e );
  void spkm3vxt_( float* vel, float* dx, int* nxv, int* np, float* afit, float* ee, float* x, float* bfit, float* xx );
}

namespace mod_spk_dvsmig3d {

  struct VariableStruct {
    int   ntIn;
    int   ntUse;
    float sr;
    float apx;
    int   irfc;
    int   itzr;
    int   nnz;
    int   nz;
    float ddz;
    float dx, dy;
    float frq1, frq2;
    float oty;
    int   nxParam;              // traces per y-line (0: from line header)
    std::string lineHeader;     // header that identifies the y-line
    bool  outputVelocity;
    // Dense velocity functions: afit(1..3) = nxv, nyv, np; then per (iy,ix): iy, ix, pairs (v,z)
    std::vector<float> afit;
    int   nxv, nyv, np;
    csTraceGather* gather;
  };

  //----------------------------------------------------------------------
  // Minimal SEISPAK reader for the velocity file (same format as INPUT_SEISPAK)
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
using namespace mod_spk_dvsmig3d;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_spk_dvsmig3d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  std::string text;
  // Defaults as in edvsmig.f
  vars->apx  = 4.0f;
  vars->irfc = 0;
  vars->itzr = 5;
  vars->frq1 = 5.0f;
  vars->frq2 = 60.0f;
  vars->oty  = 0.0f;
  vars->nnz  = 0;
  vars->ddz  = 0.0f;
  vars->dx   = 0.0f;
  vars->dy   = 0.0f;
  vars->nxParam = 0;
  vars->lineHeader = "";
  vars->outputVelocity = false;
  vars->nxv = vars->nyv = vars->np = 0;
  vars->gather = NULL;

  if( param->exists("apx") )  param->getFloat( "apx", &vars->apx );
  if( param->exists("rfc") ) {
    param->getString( "rfc", &text );
    if( !text.compare("yes") || !text.compare("YES") ) vars->irfc = 1;
    else if( !text.compare("no") || !text.compare("NO") ) vars->irfc = 0;
    else writer->error("Unknown option for 'rfc': %s", text.c_str());
  }
  if( param->exists("itzr") ) param->getInt( "itzr", &vars->itzr );
  if( param->exists("frq1") ) param->getFloat( "frq1", &vars->frq1 );
  if( param->exists("frq2") ) param->getFloat( "frq2", &vars->frq2 );
  if( param->exists("dx") )   param->getFloat( "dx", &vars->dx );
  if( param->exists("dy") )   param->getFloat( "dy", &vars->dy );
  if( param->exists("ddz") )  param->getFloat( "ddz", &vars->ddz );
  if( param->exists("nnz") )  param->getInt( "nnz", &vars->nnz );
  if( param->exists("oty") )  param->getFloat( "oty", &vars->oty );
  if( param->exists("nx") )   param->getInt( "nx", &vars->nxParam );
  if( param->exists("line_header") ) param->getString( "line_header", &vars->lineHeader );
  if( param->exists("output") ) {
    param->getString( "output", &text );
    if( !text.compare("velocity") ) vars->outputVelocity = true;
    else if( text.compare("migration") ) writer->error("Unknown option for 'output': %s", text.c_str());
  }

  float a = vars->apx;
  bool apxOK = ( a == 4.0f || a == 5.0f || a == 6.0f || a == 7.0f || a >= 11.0f );
  if( !apxOK ) writer->error("APX = %g not supported in 3D. Use 4/5 (45-degree x-y-w, split), 6/7 (phase shift + split-step) or >= 11 (multi-velocity PSPI, nvels = APX-10)", a);
  if( a >= 11.0f && (int)lrintf(a-10.0f) < 2 ) writer->error("APX >= 11 requires at least 2 velocities (APX >= 12)");
  if( vars->dx <= 0.0f ) writer->error("Parameter 'dx' (distance between depth points in x) is required");
  if( vars->dy <= 0.0f ) { vars->dy = vars->dx; writer->warning("DY not given: DY = DX = %g", vars->dx); }
  if( vars->itzr < 1 ) writer->error("ITZR must be >= 1");
  if( vars->frq2 <= vars->frq1 ) writer->error("FRQ2 must be larger than FRQ1");
  if( vars->nxParam < 0 ) writer->error("NX must be > 0");

  if( vars->nxParam == 0 ) {
    if( vars->lineHeader.empty() ) {
      if( hdef->headerExists("row") ) vars->lineHeader = "row";
      else if( hdef->headerExists("spk_rec") ) vars->lineHeader = "spk_rec";
      else writer->error("Cannot find y-lines: no header 'row' or 'spk_rec'. Give 'nx' (traces per line) or 'line_header'.");
    }
    else if( !hdef->headerExists(vars->lineHeader) ) {
      writer->error("Line header '%s' does not exist", vars->lineHeader.c_str());
    }
  }

  vars->ntIn = shdr->numSamples;
  vars->sr   = shdr->sampleInt;
  vars->ntUse = vars->ntIn;
  if( param->exists("nt") ) {
    int n = 0; param->getInt( "nt", &n );
    if( n > 0 && n < vars->ntIn ) vars->ntUse = n;
  }
  if( vars->irfc == 1 ) {
    if( vars->ddz <= 0.0f || vars->nnz <= 0 ) writer->error("Depth migration (RFC YES) requires DDZ and NNZ");
  }
  else {
    if( vars->nnz <= 0 ) vars->nnz = vars->ntUse;
  }
  vars->nz = (vars->nnz - 1)/vars->itzr + 1;
  float fnyq = 500.0f / vars->sr;
  if( vars->frq2 > fnyq ) writer->warning("FRQ2 (%f) above Nyquist (%f Hz): limited to Nyquist", vars->frq2, fnyq);

  //---------------------------------------------
  // Velocity model -> dense functions (as SPK_DVSSYN3D)
  bool hasVxt = param->exists("vxt_file"), hasVel = param->exists("vel_file");
  if( hasVxt == hasVel ) writer->error("Specify one velocity source: 'vxt_file' (SEISPAK VXT list) or 'vel_file' (SEISPAK ONL v(z) file)");
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
      spkm3vxt_( &list[0], &ddd, &gnx[g], &gnp[g], &gafit[g][0], &ee[0], &x[0], &bfit[0], &xx[0] );
      if( gnx[g] <= 0 || gnp[g] <= 0 || gnx[g] > 10000 ) writer->error("VXT file: could not build velocity functions for y-line %g", ys[g]);
      if( gnx[g] != gnx[0] ) writer->error("VXT file: y-line %g has %d DPN, y-line %g has %d. All y-lines must span the same DPN range.", ys[g], gnx[g], ys[0], gnx[0]);
    }
    int nxv = gnx[0], np = 0;
    for( int g = 0; g < ngroups; g++ ) np = std::max( np, gnp[g] );
    // y-lines of the model: 1 .. last y of the VXT (linear interpolation between VXT lines)
    int nyv = std::max( 1, (int)lrintf( ys.back() ) );
    int nw = 2*np+3;
    vars->afit.assign( 3 + (size_t)nxv*nyv*nw + 1000, 0.0f );
    for( int iy = 1; iy <= nyv; iy++ ) {
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
    vars->nxv = nxv; vars->nyv = nyv; vars->np = np;
    writer->line("  VXT file: %s  (code %g, %d y-lines from y=%g to y=%g, up to %d horizons, %d DPN)",
                 fname.c_str(), v[0], ngroups, ys.front(), ys.back(), np, nxv);
  }
  else {
    std::string fname;
    param->getString( "vel_file", &fname );
    int hmt, hmr, ntv; float dzv;
    std::vector<float> vdat;
    std::string err = readSeispak( fname, &hmt, &hmr, &ntv, &dzv, &vdat );
    if( !err.empty() ) writer->error("Velocity file '%s': %s", fname.c_str(), err.c_str());
    // As ONLVELS: record = y-line, trace = DPN
    int nxv = hmt, nyv = hmr, np = ntv, nw = 2*np+3;
    if( hmt == 1 && hmr > 1 ) writer->warning("Velocity file has one trace per record (2D layout): each record is taken as one y-line with a single DPN");
    vars->afit.assign( 3 + (size_t)nxv*nyv*nw + 1000, 0.0f );
    for( int iy = 1; iy <= nyv; iy++ ) {
      for( int ix = 1; ix <= nxv; ix++ ) {
        float* xx = &vars->afit[3 + ((size_t)(iy-1)*nxv + ix-1)*nw];
        float const* tr = &vdat[((size_t)(iy-1)*nxv + ix-1)*ntv];
        xx[0] = (float)iy; xx[1] = (float)ix;
        for( int j = 0; j < np; j++ ) { xx[2+2*j] = tr[j]; xx[3+2*j] = j*dzv; }
      }
    }
    vars->nxv = nxv; vars->nyv = nyv; vars->np = np;
    writer->line("  Velocity file (ONL): %s  (%d DPN x %d y-lines, %d samples, depth step %f)", fname.c_str(), nxv, nyv, np, dzv);
  }
  vars->afit[0] = (float)vars->nxv; vars->afit[1] = (float)vars->nyv; vars->afit[2] = (float)vars->np;

  //---------------------------------------------
  shdr->numSamples = vars->nnz;
  if( vars->irfc == 1 ) {
    shdr->sampleInt = vars->ddz;
    shdr->domain = DOMAIN_XD;
  }
  vars->gather = new csTraceGather( hdef );

  char const* apxText = ( a < 6.0f ) ? "45-degree x-y-w finite difference (split)" :
                        ( a < 11.0f ) ? "phase shift + split-step correction" : "multi-velocity PSPI";
  writer->line("  APX:  %g (%s)", a, apxText);
  writer->line("  RFC:  %s (%s migration)", vars->irfc ? "YES" : "NO", vars->irfc ? "depth" : "time");
  writer->line("  DX: %f  DY: %f  ITZR: %d   FRQ1/FRQ2: %f %f Hz   OTY: %f", vars->dx, vars->dy, vars->itzr, vars->frq1, vars->frq2, vars->oty);
  if( vars->irfc ) writer->line("  DDZ: %f   NNZ: %d  (max. depth %f)", vars->ddz, vars->nnz, (vars->nnz-1)*vars->ddz);
  else writer->line("  NNZ: %d time samples", vars->nnz);
  if( vars->nxParam > 0 ) writer->line("  Y-lines: %d traces per line (parameter NX)", vars->nxParam);
  else writer->line("  Y-lines: from trace header '%s'", vars->lineHeader.c_str());
  if( vars->outputVelocity ) writer->line("  OUTPUT VELOCITY: the velocity model used by the migration is output instead of the migrated volume");
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_spk_dvsmig3d_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csTraceHeaderDef const* hdef = env->headerDef;

  if( !edef->isLastCall() ) {
    // Collect the whole volume
    traceGather->moveTracesTo( 0, traceGather->numTraces(), vars->gather );
    edef->setTracesAreWaiting();
    return;
  }

  int ntr = vars->gather->numTraces();
  if( ntr == 0 ) return;

  //--------------------------------------------------------------
  // Organise traces in y-lines (x fast, y slow)
  int nx = 0, ny = 0;
  if( vars->nxParam > 0 ) {
    nx = vars->nxParam;
    if( ntr % nx != 0 ) writer->error("Number of traces (%d) is not a multiple of NX (%d)", ntr, nx);
    ny = ntr / nx;
  }
  else {
    int id = hdef->headerIndex( vars->lineHeader );
    std::vector<int> lineLen;
    std::vector<double> lineVal;
    for( int k = 0; k < ntr; k++ ) {
      double val = vars->gather->trace(k)->getTraceHeader()->doubleValue( id );
      if( k == 0 || val != lineVal.back() ) {
        for( size_t i = 0; i < lineVal.size(); i++ ) {
          if( lineVal[i] == val ) writer->error("Line header '%s' = %g appears again after other lines: sort the data by y-line and x", vars->lineHeader.c_str(), val);
        }
        lineVal.push_back( val );
        lineLen.push_back( 0 );
      }
      lineLen.back()++;
    }
    ny = (int)lineLen.size();
    nx = lineLen[0];
    for( int i = 1; i < ny; i++ ) {
      if( lineLen[i] != nx ) {
        char msg[512];
        snprintf( msg, sizeof(msg), "Y-line %d (%s = %g) has %d traces, first line has %d. All y-lines must have the same number of traces.",
                  i+1, vars->lineHeader.c_str(), lineVal[i], lineLen[i], nx );
        writer->error( "%s", msg );
      }
    }
    writer->line("  Lines from header '%s': %g ... %g", vars->lineHeader.c_str(), lineVal.front(), lineVal.back());
  }
  writer->line("  Data: %d traces = %d y-lines x %d traces (x)", ntr, ny, nx);
  if( ny == 1 ) writer->warning("Only one y-line: use SPK_DVSMIG2D for 2D data");

  int nt = vars->ntUse;
  std::vector<float> d( (size_t)ntr*nt );
  for( int k = 0; k < ntr; k++ ) {
    float const* s = vars->gather->trace(k)->getTraceSamples();
    memcpy( &d[(size_t)k*nt], s, nt*sizeof(float) );
  }

  //--------------------------------------------------------------
  // Velocity per depth step at each trace: VXZ(nz, nx*ny), as xvel3d:
  // x location = trace in line + OTY, y location = line number (limited to the model)
  int nz = vars->nz;
  float dt = 0.001f * vars->sr;
  std::vector<float> vxz( (size_t)nz*ntr );
  int ibeta = (int)lrintf( vars->oty );
  int np = vars->np, nw = 2*np+3;
  if( vars->nxv != nx ) writer->warning("Velocity model has %d DPN but the lines have %d traces. Trace k uses DPN k+OTY (limited to 1..%d).", vars->nxv, nx, vars->nxv);
  if( vars->nyv != ny ) writer->warning("Velocity model has %d y-lines but the data has %d. Line i uses velocity line i (limited to 1..%d).", vars->nyv, ny, vars->nyv);
  for( int i = 0; i < ny; i++ ) {
    int iy = std::min( i+1, vars->nyv );
    for( int k = 0; k < nx; k++ ) {
      int ix = k + 1 + ibeta;
      if( ix < 1 ) ix = 1;
      if( ix > vars->nxv ) ix = vars->nxv;
      float* xx = &vars->afit[3 + ((size_t)(iy-1)*vars->nxv + ix-1)*nw];
      spkvcol_( xx, &np, &nz, &vars->itzr, &vars->ddz, &dt, &vars->irfc, &vxz[((size_t)i*nx + k)*nz] );
    }
  }
  float vmin = 1e30f, vmax = 0.0f;
  for( size_t i = 0; i < vxz.size(); i++ ) {
    if( !(vxz[i] > 0.0f) ) writer->error("Invalid (zero or negative) velocity after building the velocity model. Check VXT / vel_file.");
    if( vxz[i] < vmin ) vmin = vxz[i];
    if( vxz[i] > vmax ) vmax = vxz[i];
  }
  writer->line("  Velocity range in migration: %f - %f", vmin, vmax);

  //--------------------------------------------------------------
  int nnz = vars->nnz;
  std::vector<float> out( (size_t)nnz*ntr, 0.0f );
  if( vars->outputVelocity ) {
    for( int k = 0; k < ntr; k++ ) {
      for( int m = 0; m < nnz; m++ ) {
        int step = m / vars->itzr;
        if( step >= nz ) step = nz-1;
        out[(size_t)k*nnz + m] = vxz[(size_t)k*nz + step];
      }
    }
    writer->line("  Output: velocity model (no migration)");
  }
  else {
    int info[7] = {0,0,0,0,0,0,0};
    int ierr = 0;
    float frq2 = vars->frq2;
    writer->line("  Migrating %d x %d traces...", nx, ny);
    spkmig3d_( &d[0], &nx, &ny, &nt, &vars->sr, &vxz[0], &nz, &vars->itzr, &nnz, &vars->ddz, &vars->dx, &vars->dy,
               &vars->apx, &vars->irfc, &vars->frq1, &frq2, &out[0], info, &ierr );
    writer->line("  nx, ny (FFT) %d %d  nt1 %d  frequency samples %d-%d  depth steps %d  work memory %d MB",
                 info[0], info[1], info[2], info[3], info[4], info[5], info[6]);
    if( ierr == 1 ) writer->error("No frequencies to migrate between FRQ1 and FRQ2");
    if( ierr == 2 ) writer->error("Not enough memory for migration (%d MB)", info[6]);
  }

  //--------------------------------------------------------------
  traceGather->createTraces( 0, ntr, hdef, nnz );
  for( int k = 0; k < ntr; k++ ) {
    csTrace* trc = traceGather->trace(k);
    trc->getTraceHeader()->copyFrom( vars->gather->trace(k)->getTraceHeader() );
    memcpy( trc->getTraceSamples(), &out[(size_t)k*nnz], nnz*sizeof(float) );
  }
  vars->gather->freeAllTraces();
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_spk_dvsmig3d_( csParamDef* pdef ) {
  pdef->setModule( "SPK_DVSMIG3D", "3D post-stack wave-equation migration (SEISPAK DVSMIG)",
                   "Migrates the whole 3D stacked volume. Input traces must be sorted by y-line and, inside each line, by x (DPN), "
                   "all lines with the same number of traces, spaced DX in x and DY in y. "
                   "Time migration (RFC NO) or depth migration (RFC YES, output sampled at DDZ). "
                   "Parameter names and defaults follow SEISPAK edvsmig/pdvsmig. Computation is done by the original DVSMIG Fortran routines (path NY > 1)." );

  pdef->addParam( "apx", "Migration algorithm", NUM_VALUES_FIXED );
  pdef->addValue( "4", VALTYPE_NUMBER, "4/5: 45-degree x-y-w finite difference, split in x and y; 6/7: phase shift with mean velocity plus split-step correction; "
                  ">=11: multi-velocity PSPI with APX-10 reference velocities" );

  pdef->addParam( "rfc", "Time or depth migration", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Time migration" );
  pdef->addOption( "yes", "Depth migration (requires DDZ and NNZ)" );

  pdef->addParam( "dx", "Distance between depth points in x (DDD)", NUM_VALUES_FIXED, "Same length unit as velocities (m or ft)" );
  pdef->addValue( "", VALTYPE_NUMBER, "DX" );
  pdef->addParam( "dy", "Distance between y-lines", NUM_VALUES_FIXED, "Default: DY = DX" );
  pdef->addValue( "", VALTYPE_NUMBER, "DY" );

  pdef->addParam( "nx", "Traces per y-line", NUM_VALUES_FIXED, "If not given, y-lines are found from the line header" );
  pdef->addValue( "", VALTYPE_NUMBER, "NX" );
  pdef->addParam( "line_header", "Trace header that identifies the y-line", NUM_VALUES_FIXED,
                  "A new line starts when the value changes. Default: 'row' (SPK_DVSSYN3D), else 'spk_rec' (INPUT_SEISPAK, 3D file: record = line)" );
  pdef->addValue( "row", VALTYPE_STRING, "Header name" );

  pdef->addParam( "ddz", "Depth sample rate (RFC YES)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "DDZ" );

  pdef->addParam( "nnz", "Number of output samples", NUM_VALUES_FIXED, "Depth samples for RFC YES. For time migration, default = number of input samples" );
  pdef->addValue( "0", VALTYPE_NUMBER, "NNZ" );

  pdef->addParam( "nt", "Number of input time samples to use", NUM_VALUES_FIXED, "0 = full trace" );
  pdef->addValue( "0", VALTYPE_NUMBER, "NT" );

  pdef->addParam( "itzr", "Coarseness ratio: output samples per downward-continuation step", NUM_VALUES_FIXED );
  pdef->addValue( "5", VALTYPE_NUMBER, "ITZR" );

  pdef->addParam( "frq1", "Low frequency of data", NUM_VALUES_FIXED );
  pdef->addValue( "5", VALTYPE_NUMBER, "Hz" );
  pdef->addParam( "frq2", "High frequency of data", NUM_VALUES_FIXED );
  pdef->addValue( "60", VALTYPE_NUMBER, "Hz" );

  pdef->addParam( "oty", "Shift of velocity field in x (traces)", NUM_VALUES_FIXED, "Velocity for trace k of a line is taken at DPN k+OTY" );
  pdef->addValue( "0", VALTYPE_NUMBER, "OTY" );

  pdef->addParam( "vxt_file", "File with SEISPAK VXT velocity list (horizons), as in the DVSMIG jobdeck", NUM_VALUES_FIXED,
                  "Numbers after the keyword VXT (or the whole file): code, then for each y-line: y, horizons as triplets "
                  "(velocity, DPN, depth or time) each ended by a value <= 0.5, and 0 at the end of the line. "
                  "Each y-line is read by the original XVEVENTS routine; velocities are interpolated linearly in y between VXT lines. "
                  "Data line i (1 = first) uses model line y = i." );
  pdef->addValue( "", VALTYPE_STRING, "File name" );

  pdef->addParam( "vel_file", "SEISPAK file with interval velocity v(z)", NUM_VALUES_FIXED,
                  "As option VEL -3 (ONLVELS) of DVSMIG: record = y-line, trace = DPN; the sample interval of the file is the depth step" );
  pdef->addValue( "", VALTYPE_STRING, "File name" );

  pdef->addParam( "output", "What to output", NUM_VALUES_FIXED );
  pdef->addValue( "migration", VALTYPE_OPTION );
  pdef->addOption( "migration", "Migrated volume" );
  pdef->addOption( "velocity", "QC: interval velocity used in each downward-continuation step, on the output time/depth axis" );
}

//************************************************************************************************
bool start_exec_mod_spk_dvsmig3d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

void cleanup_mod_spk_dvsmig3d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( vars->gather != NULL ) { delete vars->gather; vars->gather = NULL; }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_spk_dvsmig3d_( csParamDef* pdef ) {
  params_mod_spk_dvsmig3d_( pdef );
}
extern "C" void _init_mod_spk_dvsmig3d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spk_dvsmig3d_( param, env, writer );
}
extern "C" bool _start_exec_mod_spk_dvsmig3d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spk_dvsmig3d_( env, writer );
}
extern "C" void _exec_mod_spk_dvsmig3d_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_dvsmig3d_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_dvsmig3d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spk_dvsmig3d_( env, writer );
}
