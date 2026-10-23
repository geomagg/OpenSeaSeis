/* Copyright (c) 2026 M. A. Gallotti Guimaraes (COPPE-UFRJ).      */
/* Distributed under the OpenSeaSeis license (see LICENSE).       */
/* Module SPK_WBMUTE: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_WBMUTE
 *
 * SEISPAK program wbmute (wbmute.f): water-bottom mute of a depth section.
 * For each trace the water-bottom depth WBZ is the first sample of its velocity
 * trace (file ONLVELS, SEISPAK ONL v(z)) with v > WV+10, WBZ = (k-1)*SMP. The data
 * are scaled by a = 1 - RAMP*(WBZ - z), limited to 0..1: zero above WBZ - 1/RAMP,
 * linear ramp up to WBZ, unchanged below (wbmute.f: WV 1490, RAMP 0.02 -> 50 m).
 * Trace j of record i uses velocity trace min(j,HMT) of record min(i,HMR).
 * SEASEIS: j and i come from trace headers (default trcno, record 1); a velocity
 * trace without water bottom keeps the WBZ of the previous trace (0 at the start).
 */

namespace mod_spk_wbmute {
  struct VariableStruct {
    int hmt, hmr, ntv; float dzv;
    std::vector<float> vel;
    std::vector<float> wbz;          // per velocity trace (record-major)
    float wv, ramp;
    int idTrc, idRec;
    int nt; float sr;
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

  // SEASEIS: decimate the velocity grid (extraction to the migration grid): DPN 1, 1+SX, 1+2SX ...; lines 1, 1+SY ...
  void decimateVel( std::vector<float>* vd, int* hmt, int* hmr, int ntv, int sx, int sy ) {
    if( sx <= 1 && sy <= 1 ) return;
    int nx2 = (*hmt - 1)/sx + 1, ny2 = (*hmr - 1)/sy + 1;
    std::vector<float> out( (size_t)nx2*ny2*ntv );
    for( int iy = 0; iy < ny2; iy++ ) for( int ix = 0; ix < nx2; ix++ )
      memcpy( &out[((size_t)iy*nx2 + ix)*ntv], &(*vd)[((size_t)(iy*sy)*(*hmt) + ix*sx)*ntv], ntv*sizeof(float) );
    vd->swap( out ); *hmt = nx2; *hmr = ny2;
  }
}
using namespace mod_spk_wbmute;

void init_mod_spk_wbmute_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* v = new VariableStruct();
  edef->setVariables( v );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  std::string fname, name;
  param->getString( "vel_file", &fname );
  v->wv = 1490.0f; v->ramp = 0.02f;
  if( param->exists("wv") )   param->getFloat( "wv", &v->wv );
  if( param->exists("ramp") ) param->getFloat( "ramp", &v->ramp );
  if( v->ramp <= 0.0f ) writer->error("RAMP must be > 0");
  std::string err = readSeispak( fname, &v->hmt, &v->hmr, &v->ntv, &v->dzv, &v->vel );
  if( !err.empty() ) { std::string msg = "Velocity file '" + fname + "': " + err; writer->error( "%s", msg.c_str() ); }
  {
    int sx = 1, sy = 1;
    if( param->exists("vel_step_x") ) param->getInt( "vel_step_x", &sx );
    if( param->exists("vel_step_y") ) param->getInt( "vel_step_y", &sy );
    if( sx < 1 || sy < 1 ) writer->error("vel_step_x / vel_step_y must be >= 1");
    if( sx > 1 || sy > 1 ) {
      int h0 = v->hmt, r0 = v->hmr;
      decimateVel( &v->vel, &v->hmt, &v->hmr, v->ntv, sx, sy );
      writer->line("  Velocity grid decimated: every %d DPN and %d line(s): %d x %d -> %d x %d", sx, sy, h0, r0, v->hmt, v->hmr);
    }
  }

  int nv = v->hmt*v->hmr, nowb = 0;
  v->wbz.assign( nv, 0.0f );
  float prev = 0.0f;
  for( int k = 0; k < nv; k++ ) {
    float const* tr = &v->vel[(size_t)k*v->ntv];
    float w = prev; bool found = false;
    for( int j = 0; j < v->ntv; j++ ) if( tr[j] > v->wv + 10.0f ) { w = j*v->dzv; found = true; break; }
    if( !found ) nowb++;
    v->wbz[k] = w; prev = w;
  }
  float wmin = *std::min_element( v->wbz.begin(), v->wbz.end() ), wmax = *std::max_element( v->wbz.begin(), v->wbz.end() );

  name = "trcno"; if( param->exists("hdr_trc") ) param->getString( "hdr_trc", &name );
  if( !hdef->headerExists(name) ) writer->error("Trace header '%s' (trace number = velocity DPN) not found: set hdr_trc", name.c_str());
  v->idTrc = hdef->headerIndex( name );
  v->idRec = -1;
  if( param->exists("hdr_rec") ) {
    param->getString( "hdr_rec", &name );
    if( !hdef->headerExists(name) ) writer->error("Trace header '%s' (record number) not found", name.c_str());
    v->idRec = hdef->headerIndex( name );
  }
  v->nt = shdr->numSamples; v->sr = shdr->sampleInt;
  writer->line("  Velocity file %s: %d traces x %d records, %d samples, step %g", fname.c_str(), v->hmt, v->hmr, v->ntv, v->dzv);
  writer->line("  Water bottom (v > %g): %g - %g%s", v->wv + 10.0f, wmin, wmax, nowb ? "  (some traces without water bottom: previous kept)" : "");
  writer->line("  Mute: zero above WBZ - %g, linear ramp to WBZ (data sample interval %g)", 1.0f/v->ramp, v->sr);
}

void exec_mod_spk_wbmute_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer )
{
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csTrace* trace = traceGather->trace(0);
  csTraceHeader* th = trace->getTraceHeader();
  float* s = trace->getTraceSamples();
  int j = th->intValue( v->idTrc );
  int i = ( v->idRec >= 0 ) ? th->intValue( v->idRec ) : 1;
  int jj = std::max( 1, std::min( j, v->hmt ) ), ii = std::max( 1, std::min( i, v->hmr ) );
  float wbz = v->wbz[ (size_t)(ii-1)*v->hmt + (jj-1) ];
  for( int k = 0; k < v->nt; k++ ) {
    float z = k*v->sr;
    float a = 1.0f - v->ramp*( wbz - z );
    s[k] *= std::min( 1.0f, std::max( a, 0.0f ) );
  }
}

void params_mod_spk_wbmute_( csParamDef* pdef ) {
  pdef->setModule( "SPK_WBMUTE", "Water-bottom mute of a depth section (SEISPAK wbmute)",
                   "Water-bottom depth of each trace from its v(z) trace (first velocity > WV+10); "
                   "zero above WBZ - 1/RAMP and linear ramp to WBZ." );
  pdef->addParam( "vel_file", "SEISPAK ONL file with v(z), one trace per DPN (as ONLVELS)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addParam( "vel_step_x", "Velocity grid: DPN per output trace (extraction: output trace i uses DPN (i-1)*STEP+1)", NUM_VALUES_FIXED,
                  "Ex.: model at 25 m, migration grid at 50 m, same origin and azimuth: 2" );
  pdef->addValue( "1", VALTYPE_NUMBER, "Step" );
  pdef->addParam( "vel_step_y", "Velocity grid: lines per output line", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "Step" );
  pdef->addParam( "wv", "Water velocity: water bottom = first sample with v > WV+10", NUM_VALUES_FIXED );
  pdef->addValue( "1490", VALTYPE_NUMBER, "WV" );
  pdef->addParam( "ramp", "Ramp slope per depth unit (mute length above WBZ = 1/RAMP)", NUM_VALUES_FIXED );
  pdef->addValue( "0.02", VALTYPE_NUMBER, "RAMP" );
  pdef->addParam( "hdr_trc", "Trace header with the trace number j (velocity DPN)", NUM_VALUES_FIXED );
  pdef->addValue( "trcno", VALTYPE_STRING, "Header name" );
  pdef->addParam( "hdr_rec", "Trace header with the record number i (velocity record), default: record 1", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Header name" );
}

bool start_exec_mod_spk_wbmute_( csExecPhaseEnv* env, csLogWriter* writer ) { return true; }
void cleanup_mod_spk_wbmute_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  delete v; v = NULL;
}
extern "C" void _params_mod_spk_wbmute_( csParamDef* pdef ) { params_mod_spk_wbmute_( pdef ); }
extern "C" void _init_mod_spk_wbmute_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) { init_mod_spk_wbmute_( param, env, writer ); }
extern "C" bool _start_exec_mod_spk_wbmute_( csExecPhaseEnv* env, csLogWriter* writer ) { return start_exec_mod_spk_wbmute_( env, writer ); }
extern "C" void _exec_mod_spk_wbmute_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_wbmute_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_wbmute_( csExecPhaseEnv* env, csLogWriter* writer ) { cleanup_mod_spk_wbmute_( env, writer ); }
