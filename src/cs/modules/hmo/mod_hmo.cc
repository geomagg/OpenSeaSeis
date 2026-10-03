/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/*
 * Module HMO: hyperbolic moveout with a single constant velocity.
 * Author: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026.
 *
 *   t(x) = sqrt( t0^2 + x^2 / v^2 )
 *
 * Typical use in OBN data: v = water velocity. The direct arrival
 * t = sqrt(z^2 + x^2)/v becomes flat at t0 = z/v (z = source-receiver depth difference).
 *
 * type nmo    (default): time-variant correction t0 = sqrt(t^2 - x^2/v^2) (stretches the wavelet).
 * type static : each trace is shifted as a whole by dt = (sqrt(x^2+z^2) - z)/v (no stretch), so every
 *               event at a constant delay after the direct arrival (bubble, ghosts) stays aligned with it.
 */

#include "cseis_includes.h"
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <string>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

namespace mod_hmo {
  struct VariableStruct {
    double velocity;     // [m/s]
    bool   isApply;      // true: apply (flatten), false: remove (inverse)
    int    hdrId_offset;
    double stretchMax;   // max stretch (fraction), <= 0: no stretch mute
    float* buffer;
    bool   isStatic;     // true: static shift (no stretch)
    double depth;        // [m] constant depth (static), used if hdrId_depth < 0
    int    hdrId_depth;  // header with depth (static), -1: constant
    double flatTime;     // [ms] time where the event is placed (static), < 0: z/v
  };
}
using mod_hmo::VariableStruct;

//*************************************************************************************************
// Init phase
//
void init_mod_hmo_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csTraceHeaderDef* hdef = env->headerDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars   = new VariableStruct();
  edef->setVariables( vars );

  vars->velocity     = 0;
  vars->isApply      = true;
  vars->hdrId_offset = -1;
  vars->stretchMax   = 0;
  vars->buffer       = NULL;
  vars->isStatic     = false;
  vars->depth        = 0;
  vars->hdrId_depth  = -1;
  vars->flatTime     = -1;

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  if( shdr->domain != DOMAIN_XT ) {
    writer->error("HMO requires time-domain (XT) input.");
  }

  //--- velocity
  float velocity = 0;
  param->getFloat( "velocity", &velocity );
  if( velocity <= 0 ) {
    writer->error("Velocity must be positive: %f", velocity);
  }
  vars->velocity = velocity;

  //--- mode
  if( param->exists("mode") ) {
    std::string text;
    param->getString( "mode", &text );
    text = toLowerCase( text );
    if( !text.compare("apply") ) {
      vars->isApply = true;
    }
    else if( !text.compare("remove") ) {
      vars->isApply = false;
    }
    else {
      writer->error("Unknown option for parameter 'mode': '%s'", text.c_str());
    }
  }

  //--- offset header
  std::string hdrName = "offset";
  if( param->exists("hdr_offset") ) {
    param->getString( "hdr_offset", &hdrName );
  }
  if( !hdef->headerExists( hdrName ) ) {
    writer->error("Trace header '%s' does not exist. Compute it first (e.g. with $HDR_MATH or $POSCALC).", hdrName.c_str());
  }
  type_t type = hdef->headerType( hdrName );
  if( type != TYPE_FLOAT && type != TYPE_DOUBLE && type != TYPE_INT && type != TYPE_INT64 ) {
    writer->error("Trace header '%s' must be a number.", hdrName.c_str());
  }
  vars->hdrId_offset = hdef->headerIndex( hdrName );

  //--- type: nmo (time-variant) or static (shift)
  if( param->exists("type") ) {
    std::string text;
    param->getString( "type", &text );
    text = toLowerCase( text );
    if( !text.compare("nmo") ) vars->isStatic = false;
    else if( !text.compare("static") ) vars->isStatic = true;
    else writer->error("Unknown option for parameter 'type': '%s'", text.c_str());
  }
  std::string depthText = "rec_wdep";
  if( vars->isStatic ) {
    if( param->exists("depth") ) param->getString( "depth", &depthText );
    char* endp = NULL;
    double d = strtod( depthText.c_str(), &endp );
    if( endp != NULL && endp != depthText.c_str() && *endp == '\0' ) {
      vars->depth = fabs( d );
    }
    else {
      depthText = toLowerCase( depthText );
      if( !hdef->headerExists( depthText ) ) {
        writer->error("Depth header '%s' does not exist. Give the depth in m (e.g. 'depth 2094') or an existing header (parameter 'depth').", depthText.c_str());
      }
      if( hdef->headerType( depthText ) == TYPE_STRING ) writer->error("Depth header '%s' is a string header.", depthText.c_str());
      vars->hdrId_depth = hdef->headerIndex( depthText );
    }
    if( param->exists("flat_time") ) {
      float ft = 0;
      param->getFloat( "flat_time", &ft );
      vars->flatTime = ft;
    }
    if( param->exists("stretch_mute") ) writer->warning("Parameter 'stretch_mute' is ignored with 'type static' (no stretch).");
  }
  else if( param->exists("depth") || param->exists("flat_time") ) {
    writer->warning("Parameters 'depth' and 'flat_time' are only used with 'type static'.");
  }

  //--- stretch mute
  if( param->exists("stretch_mute") ) {
    float percent = 0;
    param->getFloat( "stretch_mute", &percent );
    if( percent < 0 ) writer->error("Stretch mute must be >= 0: %f", percent);
    vars->stretchMax = percent / 100.0;
  }

  vars->buffer = new float[shdr->numSamples];

  if( vars->isStatic ) {
    char ft[64];
    if( vars->flatTime >= 0 ) snprintf( ft, sizeof(ft), "%.2f ms", vars->flatTime );
    else snprintf( ft, sizeof(ft), "z/v (zero-offset time)" );
    writer->line("HMO: type = static (shift, no stretch), velocity = %.2f m/s, mode = %s, offset header = '%s', depth = %s, event placed at %s",
                 vars->velocity, vars->isApply ? "apply" : "remove", hdrName.c_str(),
                 vars->hdrId_depth >= 0 ? ( "header " + depthText ).c_str() : ( std::to_string( vars->depth ) + " m" ).c_str(), ft );
  }
  else {
    writer->line("HMO: type = nmo (time-variant), velocity = %.2f m/s, mode = %s, offset header = '%s', stretch mute = %s",
                 vars->velocity, vars->isApply ? "apply" : "remove", hdrName.c_str(),
                 vars->stretchMax > 0 ? "on" : "off");
  }
}

//*************************************************************************************************
// Exec phase
//
void exec_mod_hmo_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;

  csTrace* trace  = traceGather->trace(0);
  float* samples  = trace->getTraceSamples();
  int nSamples    = shdr->numSamples;
  double dt       = shdr->sampleInt / 1000.0;   // [s]

  double offset   = trace->getTraceHeader()->doubleValue( vars->hdrId_offset );
  double x2v2     = ( offset * offset ) / ( vars->velocity * vars->velocity );  // [s^2]

  float* out = vars->buffer;

  if( vars->isStatic ) {
    // Static shift of the whole trace: dt = t(x) - tRef, t(x) = sqrt(x^2+z^2)/v, tRef = z/v or flat_time
    double z = ( vars->hdrId_depth >= 0 ) ? fabs( trace->getTraceHeader()->doubleValue( vars->hdrId_depth ) ) : vars->depth;
    double tEvent = sqrt( offset*offset + z*z ) / vars->velocity;                      // [s]
    double tRef   = ( vars->flatTime >= 0 ) ? vars->flatTime / 1000.0 : z / vars->velocity;
    double shift  = tEvent - tRef;                                                     // [s], > 0: event moves up
    if( !vars->isApply ) shift = -shift;
    for( int isamp = 0; isamp < nSamples; isamp++ ) {
      double pos = isamp + shift / dt;
      int i0 = (int)floor( pos );
      if( pos < 0 || i0 > nSamples - 1 ) { out[isamp] = 0.0f; continue; }
      if( i0 == nSamples - 1 ) { out[isamp] = ( pos == (double)i0 ) ? samples[i0] : 0.0f; continue; }
      double w = pos - i0;
      out[isamp] = (float)( (1.0 - w) * samples[i0] + w * samples[i0+1] );
    }
    memcpy( samples, out, nSamples * sizeof(float) );
    return;
  }

  for( int isamp = 0; isamp < nSamples; isamp++ ) {
    double tOut = isamp * dt;
    double tIn;
    if( vars->isApply ) {
      // Output at t0 takes input from t(x) = sqrt(t0^2 + x^2/v^2)
      tIn = sqrt( tOut*tOut + x2v2 );
      if( vars->stretchMax > 0 && tOut > 0 ) {
        double stretch = tIn / tOut - 1.0;
        if( stretch > vars->stretchMax ) { out[isamp] = 0.0f; continue; }
      }
      else if( vars->stretchMax > 0 ) { out[isamp] = 0.0f; continue; }
    }
    else {
      // Inverse: output at t takes input from t0 = sqrt(t^2 - x^2/v^2)
      double arg = tOut*tOut - x2v2;
      if( arg < 0 ) { out[isamp] = 0.0f; continue; }
      tIn = sqrt( arg );
    }
    // Linear interpolation
    double pos = tIn / dt;
    int i0 = (int)pos;
    if( i0 >= nSamples - 1 ) {
      out[isamp] = ( i0 == nSamples - 1 && pos == (double)i0 ) ? samples[i0] : 0.0f;
      continue;
    }
    double w = pos - i0;
    out[isamp] = (float)( (1.0 - w) * samples[i0] + w * samples[i0+1] );
  }
  memcpy( samples, out, nSamples * sizeof(float) );
}

//*************************************************************************************************
// Parameter definition (this is what 'seaseis -m hmo' prints)
//
void params_mod_hmo_( csParamDef* pdef ) {
  pdef->setModule( "HMO", "Hyperbolic moveout correction with constant velocity (e.g. water velocity for OBN)" );

  pdef->addDoc("Applies t(x) = sqrt( t0^2 + x^2/v^2 ) using a single constant velocity v and the trace offset x.");
  pdef->addDoc("mode apply : flattens events that follow this hyperbola (output at t0 = input at t(x)).");
  pdef->addDoc("mode remove: inverse operation (output at t = input at t0 = sqrt(t^2 - x^2/v^2)).");
  pdef->addDoc("OBN: with v = water velocity, the direct arrival t = sqrt(z^2 + x^2)/v is flattened at t0 = z/v.");
  pdef->addDoc("type static: instead of the time-variant correction, each trace is shifted as a whole by dt = (sqrt(x^2+z^2) - z)/v");
  pdef->addDoc("  (no stretch). The direct arrival is flattened at z/v (or at 'flat_time') and events at a constant delay after it");
  pdef->addDoc("  (source bubble, ghosts) stay aligned with it. z = node depth, a number or a trace header ('depth').");
  pdef->addDoc("Linear interpolation between samples. Time of first sample is assumed to be 0 ms.");

  pdef->addParam( "velocity", "Moveout velocity [m/s]", NUM_VALUES_FIXED, "For OBN, use the water velocity (typically 1480-1540 m/s)" );
  pdef->addValue( "", VALTYPE_NUMBER, "Velocity [m/s]" );

  pdef->addParam( "mode", "Apply or remove the hyperbolic moveout", NUM_VALUES_FIXED );
  pdef->addValue( "apply", VALTYPE_OPTION );
  pdef->addOption( "apply", "Apply HMO (flatten hyperbola)" );
  pdef->addOption( "remove", "Remove HMO (inverse correction)" );

  pdef->addParam( "type", "Type of correction", NUM_VALUES_FIXED );
  pdef->addValue( "nmo", VALTYPE_OPTION );
  pdef->addOption( "nmo", "Time-variant: t0 = sqrt(t^2 - x^2/v^2) (wavelet is stretched)" );
  pdef->addOption( "static", "Static shift per trace: dt = (sqrt(x^2+z^2) - z)/v (no stretch)" );

  pdef->addParam( "depth", "Depth z [m] for 'type static': a number, or a trace header name", NUM_VALUES_FIXED, "e.g. 2094, or rec_wdep (written by INPUT_HDF5). Absolute value is used" );
  pdef->addValue( "rec_wdep", VALTYPE_STRING, "Depth [m] or trace header name" );

  pdef->addParam( "flat_time", "Time [ms] where the corrected event is placed ('type static')", NUM_VALUES_FIXED, "Default: z/v (zero-offset time of the direct arrival)" );
  pdef->addValue( "", VALTYPE_NUMBER, "Time [ms]" );

  pdef->addParam( "hdr_offset", "Trace header containing source-receiver offset [m]", NUM_VALUES_FIXED );
  pdef->addValue( "offset", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "stretch_mute", "Stretch mute [%]", NUM_VALUES_FIXED, "Zero samples where the moveout stretch exceeds this percentage (mode apply only). Not applied if omitted" );
  pdef->addValue( "", VALTYPE_NUMBER, "Maximum stretch [%]" );
}

//************************************************************************************************
// Start exec phase
//
bool start_exec_mod_hmo_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

//************************************************************************************************
// Cleanup phase
//
void cleanup_mod_hmo_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( vars->buffer != NULL ) {
    delete [] vars->buffer;
    vars->buffer = NULL;
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_hmo_( csParamDef* pdef ) {
  params_mod_hmo_( pdef );
}
extern "C" void _init_mod_hmo_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_hmo_( param, env, writer );
}
extern "C" bool _start_exec_mod_hmo_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_hmo_( env, writer );
}
extern "C" void _exec_mod_hmo_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_hmo_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_hmo_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_hmo_( env, writer );
}
