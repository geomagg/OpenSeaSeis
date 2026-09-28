/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module SPK_DVSMIG2D: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <sys/types.h>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_DVSMIG2D
 *
 * 2D post-stack wave-equation migration from SEISPAK program DVSMIG
 * (phase shift, 45-degree finite difference, PSPI, multi-velocity PSPI),
 * in time (RFC NO) or depth (RFC YES). Parameters follow edvsmig.f / pdvsmig.
 * The whole section is collected, migrated by the Fortran code in fortran/,
 * and output with the headers of the input traces.
 */

extern "C" {
  void spkmig2d_( float* d, int* nxsave, int* nt, float* sr, float* vxz, int* nz, int* itzr, int* nnz,
                  float* ddz, float* dx, float* apx, int* irfc, float* frq1, float* frq2,
                  float* out, int* info, int* ierr );
  void spkvcol_( float* xx, int* np, int* nz, int* itzr, float* ddz, float* dt, int* irfc, float* e );
  void spkvxt_( float* vel, float* dx, int* nxv, int* nyv, int* np, float* afit, float* ee, float* x, float* bfit, float* xx );
}

namespace mod_spk_dvsmig2d {

  // Velocity function at one location: pairs (interval velocity, depth of layer bottom)
  struct VelFunc {
    int trace;                  // trace number (1..), location of the function
    std::vector<float> xx;      // xx(1),xx(2) position, then pairs, as in xvel3d
    int np;
  };

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
    float dx;
    float frq1, frq2;
    float oty;
    std::vector<VelFunc> vxt;   // velocity functions
    bool  velFromFile;
    std::string velFile;
    bool  outputVelocity;       // QC: output velocity model instead of migration
    bool  velFromVxtFile;       // VXT horizon list (SEISPAK jobdeck format)
    std::vector<float> vxtAfit; // dense velocity functions from VEVENTS
    int   vxtNx, vxtNp;
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
  // Reads all traces of a SEISPAK file. Returns error text or "".
  std::string readSeispak( std::string const& fname, int* ntrOut, int* ntOut, float* smpOut,
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
    int ntr = (int)(hmt+0.5f) * (int)(hmr+0.5f);
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
    *ntrOut = ntr; *ntOut = nt; *smpOut = smp;
    return "";
  }
}
using namespace mod_spk_dvsmig2d;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_spk_dvsmig2d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
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
  vars->velFromFile = false;
  vars->outputVelocity = false;
  vars->velFromVxtFile = false;
  vars->vxtNx = 0; vars->vxtNp = 0;
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
  if( param->exists("ddz") )  param->getFloat( "ddz", &vars->ddz );
  if( param->exists("nnz") )  param->getInt( "nnz", &vars->nnz );
  if( param->exists("oty") )  param->getFloat( "oty", &vars->oty );
  if( param->exists("output") ) {
    param->getString( "output", &text );
    if( !text.compare("velocity") ) vars->outputVelocity = true;
    else if( text.compare("migration") ) writer->error("Unknown option for 'output': %s", text.c_str());
  }

  float a = vars->apx;
  bool apxOK = ( a == 1.0f || a == 4.0f || a == 5.0f || a == 6.0f || a == 7.0f || a >= 11.0f );
  if( !apxOK ) writer->error("APX = %f not supported in 2D. Use 1 (phase shift), 4/5 (45-degree x-w), 6/7 (PSPI) or >= 11 (multi-velocity PSPI, nvels = APX-10)", a);
  if( a >= 11.0f && (int)lrintf(a-10.0f) < 2 ) writer->error("APX >= 11 requires at least 2 velocities (APX >= 12)");
  if( vars->dx <= 0.0f ) writer->error("Parameter 'dx' (distance between depth points, DDD) is required");
  if( vars->itzr < 1 ) writer->error("ITZR must be >= 1");
  if( vars->frq2 <= vars->frq1 ) writer->error("FRQ2 must be larger than FRQ1");

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
  // Velocities
  int nLines = param->getNumLines( "vxt" );
  int nSources = ( nLines > 0 ? 1 : 0 ) + ( param->exists("vel_file") ? 1 : 0 ) + ( param->exists("vxt_file") ? 1 : 0 );
  if( nSources > 1 ) writer->error("Specify only one velocity source: 'vxt', 'vxt_file' or 'vel_file'");
  if( param->exists("vxt_file") ) {
    std::string fname;
    param->getString( "vxt_file", &fname );
    vars->velFromVxtFile = true;
    float yLine = -1.0f;
    if( param->exists("vxt_line") ) param->getFloat( "vxt_line", &yLine );
    //--- read numbers after the keyword VXT (or VEL); if no keyword, the whole file
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
    if( v.size() < 5 ) writer->error("VXT file '%s': no velocity list found", fname.c_str());
    //--- split into y-line groups: [code] [negative limits] y1 <horizons...> 0 y2 <horizons...> 0 ...
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
    int ig = 0;
    if( yLine >= 0.0f ) {
      ig = -1;
      for( size_t i = 0; i < ys.size(); i++ ) if( fabs( ys[i] - yLine ) < 0.01f ) ig = (int)i;
      if( ig < 0 ) writer->error("VXT file: line y=%g not found (%d lines, y = %g ... %g)", yLine, (int)ys.size(), ys.front(), ys.back());
    }
    //--- 2D list for VEVENTS: code, limits, 1, horizons of the selected line, 0
    std::vector<float> list;
    list.push_back( v[0] );
    list.insert( list.end(), lead.begin(), lead.end() );
    list.push_back( 1.0f );
    list.insert( list.end(), bodies[ig].begin(), bodies[ig].end() );
    list.push_back( 0.0f );
    size_t m = list.size() + 60000;
    list.resize( m, 0.0f );
    std::vector<float> ee( m ), x( 20*m + 100000 ), xx( 10001 );
    size_t nafit = 3 + (size_t)10001*203;   // up to 10000 DPN and 100 horizons
    vars->vxtAfit.assign( nafit, 0.0f );
    std::vector<float> bfit( nafit );
    int nxv = 0, nyv = 0, np = 0;
    float ddd = vars->dx;
    spkvxt_( &list[0], &ddd, &nxv, &nyv, &np, &vars->vxtAfit[0], &ee[0], &x[0], &bfit[0], &xx[0] );
    if( nxv <= 0 || np <= 0 || nxv > 10000 ) writer->error("VXT file '%s': could not build velocity functions (nx=%d, events=%d)", fname.c_str(), nxv, np);
    vars->vxtNx = nxv; vars->vxtNp = np;
    writer->line("  VXT file: %s  (%d y-lines in file, using y = %g)", fname.c_str(), (int)ys.size(), ys[ig]);
    writer->line("  VXT code %g, %d horizons, %d locations (DPN)", v[0], np, nxv);
  }
  else if( param->exists("vel_file") ) {
    param->getString( "vel_file", &vars->velFile );
    vars->velFromFile = true;
    if( nLines > 0 ) writer->error("Specify either 'vel_file' or 'vxt', not both");
  }
  else if( nLines > 0 ) {
    bool timePairs = false;
    if( param->exists("vxt_domain") ) {
      param->getString( "vxt_domain", &text );
      if( !text.compare("time") ) timePairs = true;
      else if( text.compare("depth") ) writer->error("Unknown option for 'vxt_domain': %s", text.c_str());
    }
    for( int il = 0; il < nLines; il++ ) {
      int nv = param->getNumValues( "vxt", il );
      if( nv < 3 || (nv-1) % 2 != 0 ) writer->error("'vxt' line %d: expected <trace> v1 z1 v2 z2 ... (pairs velocity, depth)", il+1);
      float trc; param->getFloatAtLine( "vxt", &trc, il, 0 );
      VelFunc f;
      f.trace = (int)lrintf(trc);
      f.np = (nv-1)/2;
      f.xx.assign( 2*f.np + 4, 0.0f );
      f.xx[0] = (float)f.trace; f.xx[1] = 1.0f;
      float zPrev = 0.0f, tPrev = 0.0f;
      for( int ip = 0; ip < f.np; ip++ ) {
        float v, z;
        param->getFloatAtLine( "vxt", &v, il, 1+2*ip );
        param->getFloatAtLine( "vxt", &z, il, 2+2*ip );
        if( v <= 0.0f ) writer->error("'vxt' line %d: velocities must be > 0", il+1);
        if( timePairs ) {          // (v, t[ms]) -> (v, z): z grows by v*dt/2 in each layer
          float t = z;
          if( t <= tPrev && ip > 0 ) writer->error("'vxt' line %d: times must increase", il+1);
          z = zPrev + v*(t - tPrev)*0.0005f;
          tPrev = t;
        }
        else if( z <= zPrev && ip > 0 ) writer->error("'vxt' line %d: depths must increase", il+1);
        f.xx[2+2*ip] = v;
        f.xx[3+2*ip] = z;
        zPrev = z;
      }
      if( il > 0 && f.trace <= vars->vxt.back().trace ) writer->error("'vxt' lines must be ordered by increasing trace number");
      vars->vxt.push_back( f );
    }
  }
  else {
    if( !vars->velFromVxtFile ) writer->error("Velocity required: use 'vxt' (pairs velocity-depth per trace), 'vxt_file' (SEISPAK VXT list) or 'vel_file' (SEISPAK file with v(z) traces)");
  }

  //---------------------------------------------
  shdr->numSamples = vars->nnz;
  if( vars->irfc == 1 ) {
    shdr->sampleInt = vars->ddz;
    shdr->domain = DOMAIN_XD;
  }
  vars->gather = new csTraceGather( hdef );

  char const* apxText = ( a == 1.0f ) ? "phase shift" : ( a < 6.0f ) ? "45-degree x-w finite difference" :
                        ( a < 11.0f ) ? "PSPI" : "multi-velocity PSPI";
  writer->line("  APX:  %g (%s)", a, apxText);
  writer->line("  RFC:  %s (%s migration)", vars->irfc ? "YES" : "NO", vars->irfc ? "depth" : "time");
  writer->line("  DX: %f   ITZR: %d   FRQ1/FRQ2: %f %f Hz   OTY: %f", vars->dx, vars->itzr, vars->frq1, vars->frq2, vars->oty);
  if( vars->irfc ) writer->line("  DDZ: %f   NNZ: %d  (max. depth %f)", vars->ddz, vars->nnz, (vars->nnz-1)*vars->ddz);
  else writer->line("  NNZ: %d time samples", vars->nnz);
  if( vars->velFromVxtFile ) {}
  else if( vars->velFromFile ) writer->line("  Velocity file: %s", vars->velFile.c_str());
  else writer->line("  Velocity functions (VXT): %d", (int)vars->vxt.size());
  if( vars->outputVelocity ) writer->line("  OUTPUT VELOCITY: the velocity model used by the migration is output instead of the migrated section");
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_spk_dvsmig2d_(
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
    // Collect the whole section
    traceGather->moveTracesTo( 0, traceGather->numTraces(), vars->gather );
    edef->setTracesAreWaiting();
    return;
  }

  int nx = vars->gather->numTraces();
  if( nx == 0 ) return;
  int nt = vars->ntUse;
  std::vector<float> d( (size_t)nx*nt );
  for( int k = 0; k < nx; k++ ) {
    float const* s = vars->gather->trace(k)->getTraceSamples();
    memcpy( &d[(size_t)k*nt], s, nt*sizeof(float) );
  }

  //--------------------------------------------------------------
  // Velocity per depth step for each trace: VXZ(nz, nx)
  int nz = vars->nz;
  float dt = 0.001f * vars->sr;
  std::vector<float> vxz( (size_t)nz*nx );
  std::vector<float> e( nz );
  int ibeta = (int)lrintf( vars->oty );   // OTY: ix = k + beta (xvel3d)
  if( vars->velFromVxtFile ) {
    int np = vars->vxtNp, nxv = vars->vxtNx, nw = 2*np+3;
    writer->line("  Data: %d traces. Data trace k uses VXT location (DPN) k%+d", nx, ibeta);
    if( nxv != nx ) writer->warning("VXT has %d locations (DPN) but the data has %d traces. Velocities are linked by trace order (DPN = k+OTY, limited to 1..%d).", nxv, nx, nxv);
    for( int k = 0; k < nx; k++ ) {
      int ix = k + 1 + ibeta;
      if( ix < 1 ) ix = 1;
      if( ix > nxv ) ix = nxv;
      float* xx = &vars->vxtAfit[3 + (size_t)(ix-1)*nw];
      spkvcol_( xx, &np, &nz, &vars->itzr, &vars->ddz, &dt, &vars->irfc, &vxz[(size_t)k*nz] );
    }
  }
  else if( vars->velFromFile ) {
    int ntrV, ntV; float dzv;
    std::vector<float> vdat;
    std::string err = readSeispak( vars->velFile, &ntrV, &ntV, &dzv, &vdat );
    if( !err.empty() ) writer->error("Velocity file '%s': %s", vars->velFile.c_str(), err.c_str());
    writer->line("  Velocity file: %d traces, %d samples, depth step %f (max. depth %f)", ntrV, ntV, dzv, (ntV-1)*dzv);
    writer->line("  Data: %d traces. Data trace k uses velocity trace k%+d", nx, ibeta);
    if( ntrV != nx ) {
      writer->warning("Velocity file has %d traces but the data has %d traces. Velocities are linked by trace order: "
                      "data trace k uses velocity trace k+OTY (limited to 1..%d).", ntrV, nx, ntrV);
    }
    float zmaxMig = vars->irfc ? (vars->nnz-1)*vars->ddz : -1.0f;
    if( zmaxMig > (ntV-1)*dzv ) {
      writer->warning("Migration goes to depth %f but the velocity file only to %f: the last velocity is extended below", zmaxMig, (ntV-1)*dzv);
    }
    std::vector<float> xx( 2*ntV + 4 );
    for( int k = 0; k < nx; k++ ) {
      int ix = k + 1 + ibeta;
      if( ix < 1 ) ix = 1;
      if( ix > ntrV ) ix = ntrV;
      xx[0] = (float)ix; xx[1] = 1.0f;          // as ONLVELS branch of xvel3d
      for( int j = 0; j < ntV; j++ ) {
        xx[2+2*j] = vdat[(size_t)(ix-1)*ntV + j];
        xx[3+2*j] = j*dzv;
      }
      xx[2*ntV+2] = 0.0f;
      spkvcol_( &xx[0], &ntV, &nz, &vars->itzr, &vars->ddz, &dt, &vars->irfc, &vxz[(size_t)k*nz] );
    }
  }
  else {
    // Velocity at each VXT location, then linear interpolation between locations
    int nf = (int)vars->vxt.size();
    for( int i = 0; i < nf; i++ ) {
      if( vars->vxt[i].trace < 1 || vars->vxt[i].trace > nx ) {
        writer->warning("VXT location at trace %d is outside the data (1..%d)", vars->vxt[i].trace, nx);
      }
    }
    std::vector< std::vector<float> > ef( nf, std::vector<float>(nz) );
    for( int i = 0; i < nf; i++ ) {
      VelFunc& f = vars->vxt[i];
      spkvcol_( &f.xx[0], &f.np, &nz, &vars->itzr, &vars->ddz, &dt, &vars->irfc, &ef[i][0] );
    }
    for( int k = 0; k < nx; k++ ) {
      float trc = (float)(k + 1 + ibeta);
      float* col = &vxz[(size_t)k*nz];
      if( nf == 1 || trc <= vars->vxt[0].trace ) { memcpy( col, &ef[0][0], nz*sizeof(float) ); continue; }
      if( trc >= vars->vxt[nf-1].trace ) { memcpy( col, &ef[nf-1][0], nz*sizeof(float) ); continue; }
      int i = 0;
      while( i < nf-2 && trc > vars->vxt[i+1].trace ) i++;
      float w = (trc - vars->vxt[i].trace) / (float)(vars->vxt[i+1].trace - vars->vxt[i].trace);
      for( int n = 0; n < nz; n++ ) col[n] = (1.0f-w)*ef[i][n] + w*ef[i+1][n];
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
  std::vector<float> out( (size_t)nnz*nx, 0.0f );
  if( vars->outputVelocity ) {
    // QC: velocity of each downward-continuation step, repeated on its ITZR output samples
    for( int k = 0; k < nx; k++ ) {
      for( int m = 0; m < nnz; m++ ) {
        int step = m / vars->itzr;
        if( step >= nz ) step = nz-1;
        out[(size_t)k*nnz + m] = vxz[(size_t)k*nz + step];
      }
    }
    writer->line("  Output: velocity model (no migration)");
  }
  else {
  int info[6] = {0,0,0,0,0,0};
  int ierr = 0;
  float frq2 = vars->frq2;
  writer->line("  Migrating %d traces...", nx);
  spkmig2d_( &d[0], &nx, &nt, &vars->sr, &vxz[0], &nz, &vars->itzr, &nnz, &vars->ddz, &vars->dx,
             &vars->apx, &vars->irfc, &vars->frq1, &frq2, &out[0], info, &ierr );
  writer->line("  nx (FFT) %d  nt1 %d  frequency samples %d-%d  depth steps %d  work memory %d MB",
               info[0], info[1], info[2], info[3], info[4], info[5]);
  if( ierr == 1 ) writer->error("No frequencies to migrate between FRQ1 and FRQ2");
  if( ierr == 2 ) writer->error("Not enough memory for migration (%d MB)", info[5]);
  }

  //--------------------------------------------------------------
  // Output: migrated traces with the headers of the input traces
  traceGather->createTraces( 0, nx, hdef, nnz );
  for( int k = 0; k < nx; k++ ) {
    csTrace* trc = traceGather->trace(k);
    trc->getTraceHeader()->copyFrom( vars->gather->trace(k)->getTraceHeader() );
    memcpy( trc->getTraceSamples(), &out[(size_t)k*nnz], nnz*sizeof(float) );
  }
  vars->gather->freeAllTraces();
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_spk_dvsmig2d_( csParamDef* pdef ) {
  pdef->setModule( "SPK_DVSMIG2D", "2D post-stack wave-equation migration (SEISPAK DVSMIG)",
                   "Migrates the whole 2D stacked section (all traces of the flow, in input order, spaced DX apart). "
                   "Time migration (RFC NO) or depth migration (RFC YES, output sampled at DDZ). "
                   "Parameter names and defaults follow SEISPAK edvsmig/pdvsmig. Computation is done by the original DVSMIG Fortran routines." );

  pdef->addParam( "apx", "Migration algorithm", NUM_VALUES_FIXED );
  pdef->addValue( "4", VALTYPE_NUMBER, "1: phase shift (vertical), 4/5: 45-degree x-w, 6/7: phase shift plus interpolation (PSPI), "
                  ">=11: multi-velocity PSPI with APX-10 reference velocities" );

  pdef->addParam( "rfc", "Time or depth migration", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Time migration" );
  pdef->addOption( "yes", "Depth migration (requires DDZ and NNZ)" );

  pdef->addParam( "dx", "Distance between depth points (DDD)", NUM_VALUES_FIXED, "Same length unit as velocities (m or ft)" );
  pdef->addValue( "", VALTYPE_NUMBER, "DX" );

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

  pdef->addParam( "oty", "Shift of velocity field in traces", NUM_VALUES_FIXED, "Velocity for trace k is taken at location k+OTY" );
  pdef->addValue( "0", VALTYPE_NUMBER, "OTY" );

  pdef->addParam( "vxt", "Velocity function at one trace location (one line per location)", NUM_VALUES_VARIABLE,
                  "Pairs of interval velocity and depth of the bottom of the layer, as in VELOCITY.VXT / DVSVELS. "
                  "Velocities are interpolated linearly between locations and held constant outside." );
  pdef->addValue( "", VALTYPE_NUMBER, "Trace number (1 = first trace)" );
  pdef->addValue( "", VALTYPE_NUMBER, "v1" );
  pdef->addValue( "", VALTYPE_NUMBER, "z1 (or t1 [ms] if vxt_domain time)" );

  pdef->addParam( "vxt_domain", "Second value of each VXT pair", NUM_VALUES_FIXED );
  pdef->addValue( "depth", VALTYPE_OPTION );
  pdef->addOption( "depth", "Depth of layer bottom (as in SEISPAK)" );
  pdef->addOption( "time", "Two-way time of layer bottom [ms], converted to depth with the interval velocities" );

  pdef->addParam( "output", "What to output", NUM_VALUES_FIXED );
  pdef->addValue( "migration", VALTYPE_OPTION );
  pdef->addOption( "migration", "Migrated section" );
  pdef->addOption( "velocity", "QC: interval velocity used in each downward-continuation step, on the output time/depth axis" );

  pdef->addParam( "vxt_file", "File with SEISPAK VXT velocity list (horizons), as in the DVSMIG jobdeck", NUM_VALUES_FIXED,
                  "Numbers after the keyword VXT (or the whole file): code, then for each y-line: y, horizons as triplets "
                  "(velocity, DPN, depth or time) each ended by a value <= 0.5, and 0 at the end of the line. "
                  "Read by the original VEVENTS routine. DPN = trace number." );
  pdef->addValue( "", VALTYPE_STRING, "File name" );

  pdef->addParam( "vxt_line", "y-line of the VXT file to use (3D models)", NUM_VALUES_FIXED, "Default: first line in the file" );
  pdef->addValue( "", VALTYPE_NUMBER, "y value as written in the file" );

  pdef->addParam( "vel_file", "SEISPAK file with interval velocity v(z), one trace per location", NUM_VALUES_FIXED,
                  "As option VEL -3 (ONLVELS) of DVSMIG: the sample interval of the file is the depth step" );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
}

//************************************************************************************************
bool start_exec_mod_spk_dvsmig2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

void cleanup_mod_spk_dvsmig2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( vars->gather != NULL ) { delete vars->gather; vars->gather = NULL; }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_spk_dvsmig2d_( csParamDef* pdef ) {
  params_mod_spk_dvsmig2d_( pdef );
}
extern "C" void _init_mod_spk_dvsmig2d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spk_dvsmig2d_( param, env, writer );
}
extern "C" bool _start_exec_mod_spk_dvsmig2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spk_dvsmig2d_( env, writer );
}
extern "C" void _exec_mod_spk_dvsmig2d_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_dvsmig2d_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_dvsmig2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spk_dvsmig2d_( env, writer );
}
