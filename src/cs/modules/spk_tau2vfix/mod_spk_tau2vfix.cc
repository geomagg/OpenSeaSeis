/* Copyright (c) 2026 M. A. Gallotti Guimaraes (COPPE-UFRJ).      */
/* Distributed under the OpenSeaSeis license (see LICENSE).       */
/* Module SPK_TAU2VFIX: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_TAU2VFIX
 *
 * SEISPAK program tau2vfix: maps the tau-deviation panels of the depth focusing
 * analysis (OFFLINE4 of xvspmig1/xvspmig2/test3dbcpmig, or PREFIX_tau.sgy of the
 * focus_* options of SPK_XVSPMIG1/2 and SPK_3DBCPMIG) to velocity coherence panels,
 * with VC = VM*(1 - DTAU/(FACTOR*TAU)), VM = average velocity of the migration (VAVG).
 * Input: one ensemble per analysis point (use ENS_DEFINE). The original routines
 * XTAU2V and FACTM are in fortran/tau2vfix_sub.f. FACTOR < 0.1 (default 0): the
 * time-variant calibration FACTM (1.0 at 0 s, 0.9 at 2 s, 0.8 at 4 s, 0.7 at 6 s, 0.65 at 8 s, 0.6 at 10 s) is used.
 *
 * Changes to the original, marked SEASEIS:
 *  - the VAVG trace of analysis point i is used (the original always read the first
 *    VAVG trace, for every record; vtrack already used trace min(i, ntraces)).
 *  - the zero line (gam = 0) is only written inside the panel.
 */

extern "C" void xtau2v_( float* tdata, float* vdata, float* vel, float* vmin, float* vmax, int* nt, float* dt, int* jnt );

namespace mod_spk_tau2vfix {
  struct VariableStruct {
    float vmin, vmax, factor;
    int nt, ntrv, ntv;
    std::vector<float> vavg;   // ntrv traces x ntv samples
    int iens;
  };

  static int const BASE = 12;
  static inline void swap4( char* b ) { char c=b[0];b[0]=b[3];b[3]=c; c=b[1];b[1]=b[2];b[2]=c; }
  static inline float getF( char const* p, bool sw ) { char b[4]; memcpy(b,p,4); if(sw) swap4(b); float v; memcpy(&v,b,4); return v; }
  static inline int   getI( char const* p, bool sw ) { char b[4]; memcpy(b,p,4); if(sw) swap4(b); int v; memcpy(&v,b,4); return v; }
  static inline short getS( char const* p, bool sw ) { char b[2]; memcpy(b,p,2); if(sw){ char c=b[0];b[0]=b[1];b[1]=c; } short v; memcpy(&v,b,2); return v; }
  static float ibm2f( unsigned int u ) {
    if( (u & 0x7fffffff) == 0 ) return 0.0f;
    int sgn = (u >> 31) ? -1 : 1, ex = (u >> 24) & 0x7f; double m = (u & 0x00ffffff) / 16777216.0;
    return (float)( sgn * m * pow( 16.0, ex - 64 ) );
  }

  // VAVG file: SEISPAK ONL ('ICE', NBYTE 2 or 4) or SEG-Y (format 1 IBM or 5 IEEE)
  static std::string readVavg( std::string const& fname, int* ntrOut, int* ntOut, std::vector<float>* data ) {
    FILE* fp = fopen( fname.c_str(), "rb" );
    if( fp == NULL ) return "cannot open file";
    std::vector<char> raw; char buf[65536]; size_t n;
    while( (n = fread(buf,1,sizeof(buf),fp)) > 0 ) raw.insert( raw.end(), buf, buf+n );
    fclose(fp);
    if( raw.size() > 200 && strncmp(&raw[0],"ICE",3) == 0 ) {
      char const* d = &raw[BASE];
      bool sw = false;
      int lh = getI(d+4,false);
      if( !(lh >= 0 && lh < 10000) ) { sw = true; lh = getI(d+4,true); }
      int ifw0 = getI(d,sw), nt = getI(d+8,sw), nbyte = getI(d+32,sw);
      float hmt = getF(d+16,sw), hmr = getF(d+20,sw);
      if( nbyte == 0 ) nbyte = 4;
      int ntr = (int)(hmt+0.5f)*(int)(hmr+0.5f);
      if( nt <= 0 || ntr <= 0 || ifw0 < 21 ) return "implausible SEISPAK header directory";
      long long rec;
      if( nbyte == 4 ) rec = 4LL*(lh+nt);
      else if( nbyte == 2 ) rec = 4LL*(lh + nt/2 + ((nt%2==0)?1:2));
      else return "only NBYTE 2 and 4 supported";
      long long first = BASE + 4LL*(ifw0-1);
      if( first + rec*ntr > (long long)raw.size() ) return "SEISPAK file truncated";
      data->assign( (size_t)ntr*nt, 0.0f );
      for( int k = 0; k < ntr; k++ ) {
        char const* p = &raw[first + rec*k] + 4*lh;
        float* out = &(*data)[(size_t)k*nt];
        if( nbyte == 4 ) { for( int i = 0; i < nt; i++ ) out[i] = getF(p+4*i,sw); }
        else {
          float fac = getF(p,sw);
          char const* s = p + 4 + ((nt%2) ? 2 : 0);
          for( int i = 0; i < nt; i++ ) out[i] = (fac > 1e-30f) ? getS(s+2*i,sw)/fac : 0.0f;
        }
      }
      *ntrOut = ntr; *ntOut = nt;
      return "";
    }
    if( raw.size() < 3600+240 ) return "neither SEISPAK nor SEG-Y";
    int nt  = (unsigned char)raw[3220]*256 + (unsigned char)raw[3221];
    int fmt = (unsigned char)raw[3224]*256 + (unsigned char)raw[3225];
    if( nt <= 0 || (fmt != 1 && fmt != 5) ) return "SEG-Y: only formats 1 (IBM) and 5 (IEEE) supported";
    long long rec = 240 + 4LL*nt;
    int ntr = (int)( (raw.size() - 3600) / rec );
    data->assign( (size_t)ntr*nt, 0.0f );
    for( int k = 0; k < ntr; k++ ) {
      unsigned char const* p = (unsigned char const*)&raw[3600 + rec*k + 240];
      for( int i = 0; i < nt; i++ ) {
        unsigned int u = ((unsigned)p[4*i]<<24) | ((unsigned)p[4*i+1]<<16) | ((unsigned)p[4*i+2]<<8) | p[4*i+3];
        float f;
        if( fmt == 5 ) memcpy( &f, &u, 4 ); else f = ibm2f( u );
        (*data)[(size_t)k*nt+i] = f;
      }
    }
    *ntrOut = ntr; *ntOut = nt;
    return "";
  }
}
using namespace mod_spk_tau2vfix;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_spk_tau2vfix_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader*  shdr = env->superHeader;
  VariableStruct* v = new VariableStruct();
  edef->setVariables( v );
  edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );

  v->nt = shdr->numSamples;
  v->iens = 0;
  param->getFloat( "vel", &v->vmin, 0 );
  param->getFloat( "vel", &v->vmax, 1 );
  if( v->vmax <= v->vmin ) writer->error("vel: VMAX must be larger than VMIN");
  v->factor = 0.0f;
  if( param->exists("factor") ) param->getFloat( "factor", &v->factor );

  std::string fn;
  param->getString( "vavg_file", &fn );
  std::string err = readVavg( fn, &v->ntrv, &v->ntv, &v->vavg );
  if( !err.empty() ) writer->error("vavg_file '%s': %s", fn.c_str(), err.c_str());
  writer->line("  VAVG file %s: %d traces x %d samples", fn.c_str(), v->ntrv, v->ntv);
  writer->line("  Velocity range %g - %g m/s", v->vmin, v->vmax);
  if( v->factor < 0.1f ) writer->line("  Factor: time-variant calibration FACTM (1.0 at 0 s ... 0.6 at 10 s)");
  else                   writer->line("  Factor: %g", v->factor);
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_spk_tau2vfix_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer )
{
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  int ntr = traceGather->numTraces();
  if( ntr < 2 ) return;
  v->iens++;
  int nt = v->nt;
  // SEASEIS: VAVG trace of this analysis point (original: always the first)
  int iv = std::min( v->iens, v->ntrv ) - 1;
  std::vector<float> vel( nt );
  for( int k = 0; k < nt; k++ ) vel[k] = v->vavg[(size_t)iv*v->ntv + std::min( k, v->ntv-1 )];
  std::vector<float> tdata( (size_t)nt*ntr ), vdata( (size_t)nt*ntr, 0.0f );
  for( int j = 0; j < ntr; j++ ) memcpy( &tdata[(size_t)j*nt], traceGather->trace(j)->getTraceSamples(), nt*sizeof(float) );
  vdata[0] = v->factor;   // as the original: factor passed in vdata(1,1)
  float dt = 0.001f*env->superHeader->sampleInt;
  xtau2v_( &tdata[0], &vdata[0], &vel[0], &v->vmin, &v->vmax, &nt, &dt, &ntr );
  for( int j = 0; j < ntr; j++ ) memcpy( traceGather->trace(j)->getTraceSamples(), &vdata[(size_t)j*nt], nt*sizeof(float) );
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_spk_tau2vfix_( csParamDef* pdef ) {
  pdef->setModule( "SPK_TAU2VFIX", "Tau-deviation panels to velocity coherence panels (SEISPAK tau2vfix)",
                   "Depth focusing analysis: maps each tau-deviation panel (traces = tau lag, samples = two-way time) to a "
                   "velocity coherence panel (traces = velocity VMIN..VMAX) with VC = VM*(1 - DTAU/(FACTOR*TAU)). "
                   "VM is the average velocity of the migration (VAVG). One ensemble per analysis point: use ENS_DEFINE. "
                   "The output has the same number of traces as the input." );
  pdef->addParam( "vavg_file", "VAVG file: one trace per analysis point (SEISPAK ONL, or PREFIX_vavg.sgy of the focus_* options)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addParam( "vel", "Velocity range of the coherence panel", NUM_VALUES_FIXED, "As 'vmin vmax' read by tau2vfix" );
  pdef->addValue( "", VALTYPE_NUMBER, "VMIN [m/s]" );
  pdef->addValue( "", VALTYPE_NUMBER, "VMAX [m/s]" );
  pdef->addParam( "factor", "Calibration factor (DTAU = FACTOR*TAU*(1-VC/VM))", NUM_VALUES_FIXED,
                  "< 0.1: time-variant FACTM table of the original (default)" );
  pdef->addValue( "0", VALTYPE_NUMBER, "FACTOR" );
}

//************************************************************************************************
bool start_exec_mod_spk_tau2vfix_( csExecPhaseEnv* env, csLogWriter* writer ) { return true; }

void cleanup_mod_spk_tau2vfix_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( v != NULL ) writer->line("  SPK_TAU2VFIX: %d analysis points", v->iens);
  delete v;
}

extern "C" void _params_mod_spk_tau2vfix_( csParamDef* pdef ) { params_mod_spk_tau2vfix_( pdef ); }
extern "C" void _init_mod_spk_tau2vfix_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) { init_mod_spk_tau2vfix_( param, env, writer ); }
extern "C" bool _start_exec_mod_spk_tau2vfix_( csExecPhaseEnv* env, csLogWriter* writer ) { return start_exec_mod_spk_tau2vfix_( env, writer ); }
extern "C" void _exec_mod_spk_tau2vfix_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_tau2vfix_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_tau2vfix_( csExecPhaseEnv* env, csLogWriter* writer ) { cleanup_mod_spk_tau2vfix_( env, writer ); }
