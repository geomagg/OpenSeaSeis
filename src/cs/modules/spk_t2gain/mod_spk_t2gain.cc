/* Copyright (c) 2026 M. A. Gallotti Guimaraes (COPPE-UFRJ).      */
/* Distributed under the OpenSeaSeis license (see LICENSE).       */
/* Module SPK_T2GAIN: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstring>
#include <cmath>
#include <vector>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_T2GAIN
 *
 * SEISPAK program t2gain.f: every sample is multiplied by t^2, t = (k-1)*dt in
 * seconds (spherical divergence correction). The original is a main program
 * with a single loop, rewritten here in C++.
 */

namespace mod_spk_t2gain {
  struct VariableStruct {
    std::vector<float> gain;
  };
}
using mod_spk_t2gain::VariableStruct;

void init_mod_spk_t2gain_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader*  shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  int nt = shdr->numSamples;
  float dt = 0.001f * shdr->sampleInt;
  vars->gain.resize( nt );
  for( int k = 0; k < nt; k++ ) {
    float t = k*dt;
    vars->gain[k] = t*t;
  }
  writer->line("  dt = %g s: data multiplied by t**2 (gain at last sample: %g)", dt, vars->gain[nt-1]);
}

void exec_mod_spk_t2gain_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  float* s = traceGather->trace(0)->getTraceSamples();
  int nt = (int)vars->gain.size();
  for( int k = 0; k < nt; k++ ) s[k] *= vars->gain[k];
}

void params_mod_spk_t2gain_( csParamDef* pdef ) {
  pdef->setModule( "SPK_T2GAIN", "Multiply data by t**2 (SEISPAK t2gain)",
                   "Every sample is multiplied by t^2, t in seconds from the first sample (t = 0). No parameters." );
}

bool start_exec_mod_spk_t2gain_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

void cleanup_mod_spk_t2gain_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  delete vars; vars = NULL;
}

extern "C" void _params_mod_spk_t2gain_( csParamDef* pdef ) {
  params_mod_spk_t2gain_( pdef );
}
extern "C" void _init_mod_spk_t2gain_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spk_t2gain_( param, env, writer );
}
extern "C" bool _start_exec_mod_spk_t2gain_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spk_t2gain_( env, writer );
}
extern "C" void _exec_mod_spk_t2gain_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_t2gain_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_t2gain_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spk_t2gain_( env, writer );
}
