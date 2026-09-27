/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module SPK_AGC: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstring>
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_AGC
 *
 * Automatic gain control from SEISPAK program xagc.f.
 * The computation is done by the original Fortran subroutine (fortran/spkagc.f);
 * this C++ file only passes each trace to it.
 */

// Fortran subroutine: all arguments by reference
extern "C" void spkagc_( float* s, float* buf1, int* nt, float* sr,
                         float* windowl, float* valmn, float* ampmax );

namespace mod_spk_agc {
  struct VariableStruct {
    float window;   // [ms]
    float valmn;
    float ampmax;
    float* buffer;
  };
}
using mod_spk_agc::VariableStruct;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_spk_agc_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader*  shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->window = 1000.0f;
  vars->valmn  = 1.0e-6f;
  vars->ampmax = 1000.0f;
  vars->buffer = NULL;

  if( param->exists("window") ) param->getFloat( "window", &vars->window );
  if( param->exists("ampmax") ) param->getFloat( "ampmax", &vars->ampmax );
  if( param->exists("valmn") )  param->getFloat( "valmn",  &vars->valmn );

  if( vars->window <= 0.0f ) writer->error("AGC window length must be > 0 ms (given: %f)", vars->window);
  if( vars->window < 2.0f*shdr->sampleInt ) {
    writer->warning("AGC window (%f ms) is shorter than 2 samples (%f ms)", vars->window, 2.0f*shdr->sampleInt);
  }
  if( shdr->numSamples < 5 ) writer->error("SPK_AGC requires at least 5 samples per trace");

  vars->buffer = new float[shdr->numSamples];

  writer->line("  AGC window [ms]:   %f  (half window: %d samples)", vars->window, (int)( 0.5f*vars->window/shdr->sampleInt ));
  writer->line("  Output amplitude:  %f", vars->ampmax);
  writer->line("  Threshold (valmn): %g", vars->valmn);
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_spk_agc_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;

  float* samples = traceGather->trace(0)->getTraceSamples();
  int   nt = shdr->numSamples;
  float sr = shdr->sampleInt;
  spkagc_( samples, vars->buffer, &nt, &sr, &vars->window, &vars->valmn, &vars->ampmax );
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_spk_agc_( csParamDef* pdef ) {
  pdef->setModule( "SPK_AGC", "Automatic gain control (SEISPAK xagc)",
                   "Each sample is divided by the mean absolute amplitude in a sliding window centred on it, and multiplied by 'ampmax'. "
                   "Leading and trailing samples with |amplitude| <= 'valmn' are left at zero, as are the first two and last two samples of the trace. "
                   "Computation is done by the original SEISPAK Fortran routine." );

  pdef->addParam( "window", "AGC window length", NUM_VALUES_FIXED );
  pdef->addValue( "1000", VALTYPE_NUMBER, "Window length [ms]" );

  pdef->addParam( "ampmax", "Output mean amplitude", NUM_VALUES_FIXED );
  pdef->addValue( "1000", VALTYPE_NUMBER, "Mean absolute amplitude after AGC" );

  pdef->addParam( "valmn", "Amplitude threshold / stabilisation", NUM_VALUES_FIXED,
                  "Samples at start/end of trace with |amplitude| <= valmn are treated as zero (muted zone). Also added to the denominator." );
  pdef->addValue( "0.000001", VALTYPE_NUMBER, "Threshold" );
}

//************************************************************************************************
// Start exec phase
//*************************************************************************************************
bool start_exec_mod_spk_agc_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

//************************************************************************************************
// Cleanup phase
//*************************************************************************************************
void cleanup_mod_spk_agc_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( vars->buffer != NULL ) { delete [] vars->buffer; vars->buffer = NULL; }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_spk_agc_( csParamDef* pdef ) {
  params_mod_spk_agc_( pdef );
}
extern "C" void _init_mod_spk_agc_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spk_agc_( param, env, writer );
}
extern "C" bool _start_exec_mod_spk_agc_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spk_agc_( env, writer );
}
extern "C" void _exec_mod_spk_agc_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_agc_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_agc_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spk_agc_( env, writer );
}
