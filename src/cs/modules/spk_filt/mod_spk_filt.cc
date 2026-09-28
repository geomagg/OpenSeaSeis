/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module SPK_FILT: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstring>
#include <cmath>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_FILT
 *
 * Zero-phase trapezoidal band-pass filter from SEISPAK program xxfilt.f.
 * Filter design, FFT length (fac235) and FFT (FFTPACK cfftf/cfftb) are done by
 * the original Fortran code in fortran/; this file only passes each trace to it.
 */

// Fortran subroutines: all arguments by reference
extern "C" {
  void spkfnt_( int* nt, int* nt1 );
  void spkfdes_( int* nt, float* sr, float* f1, float* f2, float* f3, float* f4,
                 int* nt1, float* filt, float* wsave, int* ifrq, int* ierr );
  void spkfapp_( float* s, int* nt, int* nt1, float* filt, float* a, float* wsave );
}

namespace mod_spk_filt {
  struct VariableStruct {
    int nt1;
    float* filt;
    float* wsave;
    float* work;
  };
}
using mod_spk_filt::VariableStruct;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_spk_filt_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader*  shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->nt1   = 0;
  vars->filt  = NULL;
  vars->wsave = NULL;
  vars->work  = NULL;

  float f[4];
  if( param->getNumValues( "freq" ) != 4 ) {
    writer->error("Parameter 'freq' requires 4 values: f1 f2 f3 f4 [Hz]");
  }
  for( int i = 0; i < 4; i++ ) param->getFloat( "freq", &f[i], i );

  float fnyq = 500.0f / shdr->sampleInt;
  if( f[0] < 0.0f ) writer->error("Frequencies must be >= 0 Hz");
  if( !( f[0] <= f[1] && f[1] <= f[2] && f[2] <= f[3] ) ) {
    writer->error("Frequencies must satisfy f1 <= f2 <= f3 <= f4 (given: %f %f %f %f)", f[0], f[1], f[2], f[3]);
  }
  if( f[3] >= fnyq ) writer->error("f4 (%f Hz) must be below the Nyquist frequency (%f Hz)", f[3], fnyq);

  int nt = shdr->numSamples;
  float sr = shdr->sampleInt;
  spkfnt_( &nt, &vars->nt1 );
  vars->filt  = new float[vars->nt1];
  vars->wsave = new float[4*vars->nt1 + 15];
  vars->work  = new float[2*vars->nt1];

  int ifrq[4];
  int ierr = 0;
  spkfdes_( &nt, &sr, &f[0], &f[1], &f[2], &f[3], &vars->nt1, vars->filt, vars->wsave, ifrq, &ierr );
  if( ierr != 0 ) {
    writer->error("Filter ramps too short for FFT length %d: frequency indices %d %d %d %d. "
                  "Increase the distance f1-f2 and f3-f4, or keep f4 below Nyquist (%f Hz).",
                  vars->nt1, ifrq[0], ifrq[1], ifrq[2], ifrq[3], fnyq);
  }

  writer->line("  Trapezoid [Hz]:     %f  %f  %f  %f", f[0], f[1], f[2], f[3]);
  writer->line("  Nyquist [Hz]:       %f", fnyq);
  writer->line("  FFT length (fac235): %d   (samples: %d)", vars->nt1, nt);
  writer->line("  Frequency indices:  %d  %d  %d  %d", ifrq[0], ifrq[1], ifrq[2], ifrq[3]);
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_spk_filt_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;
  float* samples = traceGather->trace(0)->getTraceSamples();
  int nt = shdr->numSamples;
  spkfapp_( samples, &nt, &vars->nt1, vars->filt, vars->work, vars->wsave );
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_spk_filt_( csParamDef* pdef ) {
  pdef->setModule( "SPK_FILT", "Trapezoidal band-pass filter (SEISPAK xxfilt)",
                   "Zero-phase band-pass filter with linear ramps between f1-f2 and f3-f4. "
                   "Traces are zero-padded to the FFT length given by fac235 (smallest even number >= nsamples with factors 2, 3, 5). "
                   "Computation is done by the original SEISPAK Fortran routines (FFTPACK)." );

  pdef->addParam( "freq", "Trapezoid corner frequencies", NUM_VALUES_FIXED,
                  "Zero below f1, linear ramp f1-f2, pass band f2-f3, linear ramp f3-f4, zero above f4" );
  pdef->addValue( "", VALTYPE_NUMBER, "f1 [Hz]" );
  pdef->addValue( "", VALTYPE_NUMBER, "f2 [Hz]" );
  pdef->addValue( "", VALTYPE_NUMBER, "f3 [Hz]" );
  pdef->addValue( "", VALTYPE_NUMBER, "f4 [Hz]" );
}

//************************************************************************************************
// Start exec phase
//*************************************************************************************************
bool start_exec_mod_spk_filt_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

//************************************************************************************************
// Cleanup phase
//*************************************************************************************************
void cleanup_mod_spk_filt_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( vars->filt  != NULL ) { delete [] vars->filt;  vars->filt  = NULL; }
  if( vars->wsave != NULL ) { delete [] vars->wsave; vars->wsave = NULL; }
  if( vars->work  != NULL ) { delete [] vars->work;  vars->work  = NULL; }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_spk_filt_( csParamDef* pdef ) {
  params_mod_spk_filt_( pdef );
}
extern "C" void _init_mod_spk_filt_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spk_filt_( param, env, writer );
}
extern "C" bool _start_exec_mod_spk_filt_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spk_filt_( env, writer );
}
extern "C" void _exec_mod_spk_filt_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_filt_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_filt_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spk_filt_( env, writer );
}
