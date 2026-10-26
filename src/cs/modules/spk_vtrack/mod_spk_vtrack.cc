/* Copyright (c) 2026 M. A. Gallotti Guimaraes (COPPE-UFRJ).      */
/* Distributed under the OpenSeaSeis license (see LICENSE).       */
/* Module SPK_VTRACK: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

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
 * Module: SPK_VTRACK
 *
 * SEISPAK program vtrack: writes the velocity track of the migration (VAVG) onto the
 * velocity coherence panels (OFFLINE5 of the migrations, PREFIX_vcoh.sgy of the
 * focus_* options, or the output of SPK_TAU2VFIX), so the velocity used can be
 * compared with the coherence maxima. At each sample the value SCALE is written on
 * the trace of the velocity VAVG (as the original: trace iv+1, with
 * iv = 1 + nint((ntr-1)*(VAVG-VMIN)/(VMAX-VMIN))).
 * One ensemble per analysis point (use ENS_DEFINE); VAVG trace min(i, ntraces).
 *
 * Changes to the original, marked SEASEIS:
 *  - when iv = ntr the original wrote one trace beyond the panel; that sample is skipped.
 */


namespace mod_spk_vtrack {
  struct VariableStruct {
    float vmin, vmax, scale;
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
using namespace mod_spk_vtrack;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_spk_vtrack_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
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
  v->scale = -1.0f;
  if( param->exists("scale") ) param->getFloat( "scale", &v->scale );

  std::string fn;
  param->getString( "vavg_file", &fn );
  std::string err = readVavg( fn, &v->ntrv, &v->ntv, &v->vavg );
  if( !err.empty() ) writer->error("vavg_file '%s': %s", fn.c_str(), err.c_str());
  writer->line("  VAVG file %s: %d traces x %d samples", fn.c_str(), v->ntrv, v->ntv);
  writer->line("  Velocity range %g - %g m/s, track value %g", v->vmin, v->vmax, v->scale);
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_spk_vtrack_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer )
{
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  int ntr = traceGather->numTraces();
  if( ntr < 2 ) return;
  v->iens++;
  int ii = std::min( v->iens, v->ntrv ) - 1;
  for( int k = 0; k < v->nt; k++ ) {
    float v1 = v->vavg[(size_t)ii*v->ntv + std::min( k, v->ntv-1 )];
    int iv = 1 + (int)lrint( (ntr-1)*(v1 - v->vmin)/(v->vmax - v->vmin) );
    iv = std::max( 1, std::min( ntr, iv ) );
    // original: scr(k+iv*nt) = scale  -> trace iv+1 (1-based) = index iv
    if( iv < ntr ) traceGather->trace(iv)->getTraceSamples()[k] = v->scale;   // SEASEIS: iv = ntr skipped
  }
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_spk_vtrack_( csParamDef* pdef ) {
  pdef->setModule( "SPK_VTRACK", "Velocity track on velocity coherence panels (SEISPAK vtrack)",
                   "Writes SCALE on the trace of the migration average velocity (VAVG) at every sample of the velocity "
                   "coherence panels. One ensemble per analysis point: use ENS_DEFINE." );
  pdef->addParam( "vavg_file", "VAVG file: one trace per analysis point (SEISPAK ONL, or PREFIX_vavg.sgy of the focus_* options)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addParam( "vel", "Velocity range of the coherence panel", NUM_VALUES_FIXED, "As 'vmin vmax' read by vtrack" );
  pdef->addValue( "", VALTYPE_NUMBER, "VMIN [m/s]" );
  pdef->addValue( "", VALTYPE_NUMBER, "VMAX [m/s]" );
  pdef->addParam( "scale", "Value written on the track", NUM_VALUES_FIXED );
  pdef->addValue( "-1", VALTYPE_NUMBER, "SCALE" );
}

//************************************************************************************************
bool start_exec_mod_spk_vtrack_( csExecPhaseEnv* env, csLogWriter* writer ) { return true; }

void cleanup_mod_spk_vtrack_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  delete v;
}

extern "C" void _params_mod_spk_vtrack_( csParamDef* pdef ) { params_mod_spk_vtrack_( pdef ); }
extern "C" void _init_mod_spk_vtrack_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) { init_mod_spk_vtrack_( param, env, writer ); }
extern "C" bool _start_exec_mod_spk_vtrack_( csExecPhaseEnv* env, csLogWriter* writer ) { return start_exec_mod_spk_vtrack_( env, writer ); }
extern "C" void _exec_mod_spk_vtrack_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_vtrack_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_vtrack_( csExecPhaseEnv* env, csLogWriter* writer ) { cleanup_mod_spk_vtrack_( env, writer ); }
