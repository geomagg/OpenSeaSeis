/* Copyright (c) 2026 M. A. Gallotti Guimaraes (COPPE-UFRJ).      */
/* Distributed under the OpenSeaSeis license (see LICENSE).       */
/* Module SPK_LSQFIX2: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstring>
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_LSQFIX2
 *
 * SEISPAK program lsqfix2 (input: echo RNGMAX | lsqfix2), used before lsqmatch:
 *  RNGMAX > 0: keeps only traces with |RANG| < RNGMAX; in each one only the
 *    window around the direct arrival is kept:
 *    tt = 1000*sqrt( (GWDE/VEL)^2 + (RANG/VEL)^2 ) - 100 ms, data kept in [tt, tt+200]
 *    (VEL = 1500 m/s), everything else is zeroed. Output: one record.
 *  RNGMAX < 0: zeroes the traces with |RANG| < |RNGMAX| ("inside mute"), all
 *    traces are kept.
 * SEASEIS: in the RNGMAX < 0 branch the original wrote every trace at the position
 * counted only for the zeroed traces (overwriting the others); here all traces are
 * kept in order. Window start/length and VEL are parameters (defaults as the original).
 */

namespace mod_spk_lsqfix2 {
  struct VariableStruct { double rngmax, vel, wstart, wlen; int idOff, idDep; int nkept, nin; };
}
using namespace mod_spk_lsqfix2;

void init_mod_spk_lsqfix2_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* v = new VariableStruct();
  edef->setVariables( v );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  float f;
  v->rngmax = 500.; v->vel = 1500.; v->wstart = -100.; v->wlen = 200.;
  if( param->exists("rngmax") ) { param->getFloat("rngmax", &f); v->rngmax = f; }
  if( param->exists("vel") )    { param->getFloat("vel", &f);    v->vel = f; }
  if( param->exists("wstart") ) { param->getFloat("wstart", &f); v->wstart = f; }
  if( param->exists("wlen") )   { param->getFloat("wlen", &f);   v->wlen = f; }
  if( v->rngmax == 0. ) writer->error("RNGMAX must not be zero");
  string name = hdef->headerExists("spk_rang") ? "spk_rang" : "offset";
  if( param->exists("hdr_offset") ) param->getString("hdr_offset", &name);
  if( !hdef->headerExists(name) ) writer->error("Offset header '%s' not found: set hdr_offset", name.c_str());
  v->idOff = hdef->headerIndex( name );
  v->idDep = -1;
  if( v->rngmax > 0. ) {
    name = hdef->headerExists("spk_gwde") ? "spk_gwde" : "rec_wdep";
    if( param->exists("hdr_depth") ) param->getString("hdr_depth", &name);
    if( !hdef->headerExists(name) ) writer->error("Receiver water depth header '%s' (SEISPAK GWDE) not found: set hdr_depth", name.c_str());
    v->idDep = hdef->headerIndex( name );
    writer->line("  Keep |offset| < %g; window %g to %g ms around the direct arrival (VEL %g)", v->rngmax, v->wstart, v->wstart+v->wlen, v->vel);
  }
  else writer->line("  Zero traces with |offset| < %g (all traces kept)", -v->rngmax);
  v->nkept = 0; v->nin = 0;
}

void exec_mod_spk_lsqfix2_( csTraceGather* traceGather, int* port, int* numTrcToKeep,
                            csExecPhaseEnv* env, csLogWriter* writer )
{
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;
  csTrace* trace = traceGather->trace(0);
  csTraceHeader* th = trace->getTraceHeader();
  float* s = trace->getTraceSamples();
  int nt = shdr->numSamples; double sr = shdr->sampleInt;
  double rng = th->doubleValue( v->idOff );
  v->nin++;
  if( v->rngmax > 0. ) {
    if( !( fabs(rng) < v->rngmax ) ) { traceGather->freeAllTraces(); return; }
    double gzel = th->doubleValue( v->idDep );
    double t0 = gzel/v->vel;
    double tt = 1000.*sqrt( t0*t0 + pow( rng/v->vel, 2 ) ) + v->wstart;
    for( int k = 1; k <= nt; k++ ) {
      double t = (k-1)*sr;
      if( t < tt || t > tt + v->wlen ) s[k-1] = 0.0f;
    }
    v->nkept++;
  }
  else {
    if( fabs(rng) < -v->rngmax ) { for( int k = 0; k < nt; k++ ) s[k] = 0.0f; v->nkept++; }
  }
}

void params_mod_spk_lsqfix2_( csParamDef* pdef ) {
  pdef->setModule( "SPK_LSQFIX2", "Direct-arrival window of near traces (SEISPAK lsqfix2)",
                   "RNGMAX > 0: keep traces with |offset| < RNGMAX and only the window [tt-100, tt+100] ms around the direct arrival "
                   "tt = 1000*sqrt((GWDE/1500)^2+(RANG/1500)^2). RNGMAX < 0: zero the traces with |offset| < |RNGMAX|." );
  pdef->addParam( "rngmax", "RNGMAX [m] (< 0: zero the near traces instead)", NUM_VALUES_FIXED );
  pdef->addValue( "500", VALTYPE_NUMBER, "RNGMAX" );
  pdef->addParam( "vel", "Water velocity for the direct arrival [m/s]", NUM_VALUES_FIXED );
  pdef->addValue( "1500", VALTYPE_NUMBER, "VEL" );
  pdef->addParam( "wstart", "Window start relative to the direct arrival [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "-100", VALTYPE_NUMBER, "ms" );
  pdef->addParam( "wlen", "Window length [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "200", VALTYPE_NUMBER, "ms" );
  pdef->addParam( "hdr_offset", "Offset header (SEISPAK RANG)", NUM_VALUES_FIXED, "Default: spk_rang, else offset" );
  pdef->addValue( "spk_rang", VALTYPE_STRING, "Header name" );
  pdef->addParam( "hdr_depth", "Receiver water depth header (SEISPAK GWDE)", NUM_VALUES_FIXED, "Default: spk_gwde, else rec_wdep" );
  pdef->addValue( "spk_gwde", VALTYPE_STRING, "Header name" );
}
bool start_exec_mod_spk_lsqfix2_( csExecPhaseEnv* env, csLogWriter* writer ) { return true; }
void cleanup_mod_spk_lsqfix2_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( v ) writer->line("  SPK_LSQFIX2: %d of %d traces %s", v->nkept, v->nin, v->rngmax > 0 ? "kept" : "zeroed");
  delete v; v = NULL;
}
extern "C" void _params_mod_spk_lsqfix2_( csParamDef* pdef ) { params_mod_spk_lsqfix2_( pdef ); }
extern "C" void _init_mod_spk_lsqfix2_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) { init_mod_spk_lsqfix2_( param, env, writer ); }
extern "C" bool _start_exec_mod_spk_lsqfix2_( csExecPhaseEnv* env, csLogWriter* writer ) { return start_exec_mod_spk_lsqfix2_( env, writer ); }
extern "C" void _exec_mod_spk_lsqfix2_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) { exec_mod_spk_lsqfix2_( traceGather, port, numTrcToKeep, env, writer ); }
extern "C" void _cleanup_mod_spk_lsqfix2_( csExecPhaseEnv* env, csLogWriter* writer ) { cleanup_mod_spk_lsqfix2_( env, writer ); }
