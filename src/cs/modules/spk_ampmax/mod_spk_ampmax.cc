/* Copyright (c) 2026 M. A. Gallotti Guimaraes (COPPE-UFRJ).      */
/* Distributed under the OpenSeaSeis license (see LICENSE).       */
/* Module SPK_AMPMAX: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cmath>
#include <cstring>
#include <vector>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_AMPMAX
 *
 * SEISPAK program ampmax (xampmax.f): scales the whole data set so that the
 * "typical" trace maximum becomes AMPMAX. AMX1 = largest |amplitude| of all traces;
 * AVMAX = mean of the trace maxima that lie strictly between 0.90*AMX1 and
 * 0.98*AMX1; every sample is multiplied by AMPMAX/AVMAX.
 * SEASEIS: if no trace maximum lies in that band (xampmax.f divides by zero),
 * AVMAX = AMX1; an all-zero data set is left unchanged.
 */

namespace mod_spk_ampmax {
  struct VariableStruct {
    float ampmax;
    csTraceGather* gather;
  };
}
using namespace mod_spk_ampmax;

void init_mod_spk_ampmax_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* v = new VariableStruct();
  edef->setVariables( v );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  v->ampmax = 1.0f;
  if( param->exists("ampmax") ) param->getFloat( "ampmax", &v->ampmax );
  v->gather = new csTraceGather( hdef );
  writer->line("  AMPMAX %g: data scaled so that the mean of the trace maxima in (0.90, 0.98) x global max becomes AMPMAX", v->ampmax);
}

void exec_mod_spk_ampmax_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer )
{
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  if( !edef->isLastCall() ) {
    traceGather->moveTracesTo( 0, traceGather->numTraces(), v->gather );
    edef->setTracesAreWaiting();
    return;
  }
  int ntr = v->gather->numTraces();
  if( ntr == 0 ) return;
  std::vector<float> mx( ntr, 0.0f );
  float amx1 = 0.0f;
  for( int k = 0; k < ntr; k++ ) {
    csTrace* t = v->gather->trace(k);
    float const* s = t->getTraceSamples();
    int nt = t->numSamples();
    float a = 0.0f;
    for( int i = 0; i < nt; i++ ) a = std::max( a, (float)fabs( s[i] ) );
    mx[k] = a; amx1 = std::max( amx1, a );
  }
  double sum = 0.0; int kk = 0;
  for( int k = 0; k < ntr; k++ ) {
    if( 0.9f*amx1 < mx[k] && mx[k] < 0.98f*amx1 ) { sum += mx[k]; kk++; }
  }
  double avmax = ( kk > 0 ) ? sum/kk : amx1;                       // SEASEIS: no division by zero
  writer->line("  max amp found %g;  %d traces in (0.90, 0.98) x max, average max amp %g", amx1, kk, avmax);
  if( kk == 0 ) writer->warning("No trace maximum between 0.90 and 0.98 of the largest one: the largest maximum is used");
  double scale = ( avmax > 0.0 ) ? v->ampmax/avmax : 1.0;
  writer->line("  scale applied %g", scale);
  for( int k = 0; k < ntr; k++ ) {
    csTrace* t = v->gather->trace(k);
    float* s = t->getTraceSamples();
    int nt = t->numSamples();
    for( int i = 0; i < nt; i++ ) s[i] = (float)( scale*s[i] );
  }
  v->gather->moveTracesTo( 0, ntr, traceGather );
}

void params_mod_spk_ampmax_( csParamDef* pdef ) {
  pdef->setModule( "SPK_AMPMAX", "Scale the data set to a target maximum amplitude (SEISPAK ampmax)",
                   "AVMAX = mean of the trace maxima between 0.90 and 0.98 of the largest |amplitude| of all traces; "
                   "every sample is multiplied by AMPMAX/AVMAX. All traces are held in memory." );
  pdef->addParam( "ampmax", "Target amplitude (SEISPAK: ampmax FILE AMPMAX)", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "AMPMAX" );
}

bool start_exec_mod_spk_ampmax_( csExecPhaseEnv* env, csLogWriter* writer ) { return true; }
void cleanup_mod_spk_ampmax_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( v->gather != NULL ) { delete v->gather; v->gather = NULL; }
  delete v; v = NULL;
}
extern "C" void _params_mod_spk_ampmax_( csParamDef* pdef ) { params_mod_spk_ampmax_( pdef ); }
extern "C" void _init_mod_spk_ampmax_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) { init_mod_spk_ampmax_( param, env, writer ); }
extern "C" bool _start_exec_mod_spk_ampmax_( csExecPhaseEnv* env, csLogWriter* writer ) { return start_exec_mod_spk_ampmax_( env, writer ); }
extern "C" void _exec_mod_spk_ampmax_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_ampmax_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_ampmax_( csExecPhaseEnv* env, csLogWriter* writer ) { cleanup_mod_spk_ampmax_( env, writer ); }
