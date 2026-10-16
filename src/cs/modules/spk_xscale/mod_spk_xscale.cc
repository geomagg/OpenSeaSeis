/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module SPK_XSCALE: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstring>
#include <cmath>
#include <vector>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_XSCALE
 *
 * SEISPAK programs xscale.f / scale.f: multiply every trace by a constant
 * (one value, as scale.f) or by a time-variant scale interpolated linearly
 * between (amplitude, time) pairs (xscale.f). The interpolation is done by the
 * original Fortran routine XINT2 (fortran/xscale_sub.f).
 *
 * Changes to the original, marked SEASEIS:
 *  - the scale function starts at time 0 (the original started the time axis at
 *    the first pair's time, so the function was shifted when t1 > 0)
 *  - outside the pairs the first/last scale is held ('outside constant');
 *    'outside zero' keeps the original behaviour (scale 0 before the first and
 *    after the last time, i.e. those samples are zeroed)
 */

extern "C" void xint2_( float* f, float* x, float* y, float* x1, float* dx, int* n, int* m );

namespace mod_spk_xscale {
  struct VariableStruct {
    std::vector<float> scale;   // scale per sample
  };
}
using mod_spk_xscale::VariableStruct;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_spk_xscale_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader*  shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  int nt = shdr->numSamples;
  float sr = shdr->sampleInt;

  std::vector<float> v;
  int nLines = param->getNumLines( "scale" );
  for( int il = 0; il < nLines; il++ ) {
    int nv = param->getNumValues( "scale", il );
    for( int iv = 0; iv < nv; iv++ ) { float f; param->getFloatAtLine( "scale", &f, il, iv ); v.push_back( f ); }
  }
  if( v.empty() ) writer->error("Parameter 'scale' is required: one value (constant) or pairs amplitude time[ms]");

  bool outsideZero = false;
  if( param->exists("outside") ) {
    std::string text;
    param->getString( "outside", &text );
    if( !text.compare("zero") ) outsideZero = true;
    else if( text.compare("constant") ) writer->error("Unknown option for 'outside': %s", text.c_str());
  }

  vars->scale.assign( nt, 0.0f );
  if( v.size() == 1 ) {
    for( int k = 0; k < nt; k++ ) vars->scale[k] = v[0];
    writer->line("  Constant scale: %g", v[0]);
  }
  else {
    if( v.size() % 2 != 0 ) writer->error("'scale': give one value, or pairs amplitude time[ms] (got %d values)", (int)v.size());
    int nn = (int)v.size()/2;
    std::vector<float> amp( nn ), z( nn );
    for( int k = 0; k < nn; k++ ) {
      amp[k] = v[2*k];
      z[k]   = v[2*k+1];
      if( k > 0 && z[k] <= z[k-1] ) writer->error("'scale': times must increase (pair %d: %g <= %g ms)", k+1, z[k], z[k-1]);
      writer->line("  amplitude %12g  at %10.2f ms", amp[k], z[k]);
    }
    float x1 = 0.0f;   // SEASEIS: original used x1 = z(1)
    int n = nt, m = nn;
    xint2_( &vars->scale[0], &z[0], &amp[0], &x1, &sr, &n, &m );
    if( !outsideZero ) {
      // SEASEIS: hold first/last value outside the pairs (XINT2 returns 0 there)
      for( int k = 0; k < nt; k++ ) {
        float t = k*sr;
        if( t < z[0] ) vars->scale[k] = amp[0];
        else if( t > z[nn-1] ) vars->scale[k] = amp[nn-1];
      }
    }
    writer->line("  Outside the pairs: %s", outsideZero ? "scale 0 (original xscale)" : "first/last scale held");
  }
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_spk_xscale_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  float* s = traceGather->trace(0)->getTraceSamples();
  int nt = (int)vars->scale.size();
  for( int k = 0; k < nt; k++ ) s[k] *= vars->scale[k];
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_spk_xscale_( csParamDef* pdef ) {
  pdef->setModule( "SPK_XSCALE", "Constant or time-variant scaling (SEISPAK xscale / scale)",
                   "One value: every sample is multiplied by it (SEISPAK scale). Pairs amplitude time[ms]: the scale is interpolated "
                   "linearly in time (original XINT2 routine) and applied sample by sample (SEISPAK xscale)." );

  pdef->addParam( "scale", "Scale: one value, or pairs amplitude time[ms]", NUM_VALUES_VARIABLE,
                  "As the values read by xscale: a1 t1 a2 t2 ... (amplitude first). Several lines are joined." );
  pdef->addValue( "", VALTYPE_NUMBER, "Scale (or amplitude a1)" );
  pdef->addValue( "", VALTYPE_NUMBER, "Time t1 [ms]" );

  pdef->addParam( "outside", "Scale before the first and after the last pair", NUM_VALUES_FIXED );
  pdef->addValue( "constant", VALTYPE_OPTION );
  pdef->addOption( "constant", "Hold the first/last amplitude" );
  pdef->addOption( "zero", "Scale 0, samples are zeroed (as the original xscale)" );
}

//************************************************************************************************
bool start_exec_mod_spk_xscale_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

void cleanup_mod_spk_xscale_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  delete vars; vars = NULL;
}

extern "C" void _params_mod_spk_xscale_( csParamDef* pdef ) {
  params_mod_spk_xscale_( pdef );
}
extern "C" void _init_mod_spk_xscale_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spk_xscale_( param, env, writer );
}
extern "C" bool _start_exec_mod_spk_xscale_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spk_xscale_( env, writer );
}
extern "C" void _exec_mod_spk_xscale_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_xscale_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_xscale_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spk_xscale_( env, writer );
}
