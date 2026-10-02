/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/*
 * Module RMS_MAGG: RMS amplitude in a time window centred on a picked time.
 * Author: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026.
 *
 * For each trace:
 *   tp   = value of the pick header [ms] (time relative to the first sample)
 *   i1   = round( (tp - length/2) / dt ),  i2 = round( (tp + length/2) / dt )   (clipped to the trace)
 *   rms  = sqrt( sum_{i=i1..i2} (s_i - m)^2 / (i2-i1+1) ),  m = 0 or window mean (remove_mean yes)
 * The result is stored in a trace header. Samples are not changed.
 * Traces without a valid pick (pick <= 0, or window completely outside the trace) get 'null_value'.
 */

#include "cseis_includes.h"
#include <cmath>
#include <string>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

namespace mod_rms_magg {
  struct VariableStruct {
    int    hdrId_pick;
    int    hdrId_rms;
    int    hdrId_nsamp;     // -1 if not requested
    double halfLength;      // [ms]
    bool   removeMean;
    bool   zeroIsNoPick;
    float  nullValue;
    int    numTraces;
    int    numNoPick;
    int    numOutside;
    int    numClipped;
  };
}
using mod_rms_magg::VariableStruct;

//*************************************************************************************************
// Init phase
//
void init_mod_rms_magg_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csTraceHeaderDef* hdef = env->headerDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars   = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->hdrId_pick   = -1;
  vars->hdrId_rms    = -1;
  vars->hdrId_nsamp  = -1;
  vars->halfLength   = 0;
  vars->removeMean   = false;
  vars->zeroIsNoPick = true;
  vars->nullValue    = 0.0f;
  vars->numTraces = vars->numNoPick = vars->numOutside = vars->numClipped = 0;

  if( shdr->domain != DOMAIN_XT ) {
    writer->error("RMS_MAGG requires time-domain (XT) input.");
  }

  //--- Pick header
  std::string hdrPick = "time_pick";
  if( param->exists("hdr_pick") ) param->getString( "hdr_pick", &hdrPick );
  hdrPick = toLowerCase( hdrPick );
  if( !hdef->headerExists( hdrPick ) ) {
    writer->error("Pick header '%s' does not exist. Pick first (e.g. module PICKING) or set parameter 'hdr_pick'.", hdrPick.c_str());
  }
  if( hdef->headerType( hdrPick ) == TYPE_STRING ) writer->error("Pick header '%s' is a string header.", hdrPick.c_str());
  vars->hdrId_pick = hdef->headerIndex( hdrPick );

  //--- Window length (total, centred on the pick)
  float length = 0;
  param->getFloat( "length", &length );
  if( length <= 0 ) writer->error("Window length must be positive: %f ms", length);
  vars->halfLength = 0.5 * length;
  if( length < shdr->sampleInt ) {
    writer->warning("Window length (%f ms) is shorter than the sample interval (%f ms): only one sample is used.", length, shdr->sampleInt);
  }

  //--- Options
  if( param->exists("remove_mean") ) {
    std::string text;
    param->getString( "remove_mean", &text );
    text = toLowerCase( text );
    if( !text.compare("yes") ) vars->removeMean = true;
    else if( !text.compare("no") ) vars->removeMean = false;
    else writer->error("Unknown option for parameter 'remove_mean': '%s'", text.c_str());
  }
  if( param->exists("zero_pick") ) {
    std::string text;
    param->getString( "zero_pick", &text );
    text = toLowerCase( text );
    if( !text.compare("ignore") ) vars->zeroIsNoPick = true;
    else if( !text.compare("use") ) vars->zeroIsNoPick = false;
    else writer->error("Unknown option for parameter 'zero_pick': '%s'", text.c_str());
  }
  if( param->exists("null_value") ) param->getFloat( "null_value", &vars->nullValue );

  //--- Output header(s)
  std::string hdrRms = "rms_pick";
  if( param->exists("hdr_rms") ) param->getString( "hdr_rms", &hdrRms );
  hdrRms = toLowerCase( hdrRms );
  if( hdef->headerExists( hdrRms ) ) {
    if( hdef->headerType( hdrRms ) == TYPE_STRING ) writer->error("Output header '%s' exists and is a string header.", hdrRms.c_str());
  }
  else {
    hdef->addHeader( TYPE_FLOAT, hdrRms, "RMS amplitude in window centred on pick (RMS_MAGG)" );
  }
  vars->hdrId_rms = hdef->headerIndex( hdrRms );

  if( param->exists("hdr_nsamples") ) {
    std::string hdrN;
    param->getString( "hdr_nsamples", &hdrN );
    hdrN = toLowerCase( hdrN );
    if( !hdef->headerExists( hdrN ) ) hdef->addHeader( TYPE_INT, hdrN, "Number of samples used by RMS_MAGG" );
    vars->hdrId_nsamp = hdef->headerIndex( hdrN );
  }

  writer->line( "  Pick header:      %s [ms]", hdrPick.c_str() );
  writer->line( "  Window:           %g ms total (pick - %g ms ... pick + %g ms), %d samples", 2*vars->halfLength,
                vars->halfLength, vars->halfLength, (int)floor( 2*vars->halfLength / shdr->sampleInt + 0.5 ) + 1 );
  writer->line( "  Remove mean:      %s", vars->removeMean ? "yes" : "no" );
  writer->line( "  Output header:    %s", hdrRms.c_str() );
}

//*************************************************************************************************
// Exec phase
//
void exec_mod_rms_magg_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  csTrace* trace = traceGather->trace(0);
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;
  csTraceHeader* trcHdr = trace->getTraceHeader();
  float const* samples = trace->getTraceSamples();
  int nSamples = shdr->numSamples;
  double dt = shdr->sampleInt;

  vars->numTraces++;
  double pick = trcHdr->doubleValue( vars->hdrId_pick );

  if( ( vars->zeroIsNoPick && pick <= 0.0 ) || pick < 0.0 ) {
    vars->numNoPick++;
    trcHdr->setFloatValue( vars->hdrId_rms, vars->nullValue );
    if( vars->hdrId_nsamp >= 0 ) trcHdr->setIntValue( vars->hdrId_nsamp, 0 );
    return;
  }

  int i1 = (int)floor( ( pick - vars->halfLength ) / dt + 0.5 );
  int i2 = (int)floor( ( pick + vars->halfLength ) / dt + 0.5 );
  if( i2 < 0 || i1 > nSamples-1 ) {
    vars->numOutside++;
    trcHdr->setFloatValue( vars->hdrId_rms, vars->nullValue );
    if( vars->hdrId_nsamp >= 0 ) trcHdr->setIntValue( vars->hdrId_nsamp, 0 );
    return;
  }
  if( i1 < 0 || i2 > nSamples-1 ) {
    vars->numClipped++;
    if( i1 < 0 ) i1 = 0;
    if( i2 > nSamples-1 ) i2 = nSamples-1;
  }

  int n = i2 - i1 + 1;
  double mean = 0.0;
  if( vars->removeMean ) {
    for( int i = i1; i <= i2; i++ ) mean += samples[i];
    mean /= n;
  }
  double sum = 0.0;
  for( int i = i1; i <= i2; i++ ) {
    double d = samples[i] - mean;
    sum += d * d;
  }
  trcHdr->setFloatValue( vars->hdrId_rms, (float)sqrt( sum / n ) );
  if( vars->hdrId_nsamp >= 0 ) trcHdr->setIntValue( vars->hdrId_nsamp, n );
}

//*************************************************************************************************
// Parameter definition (help: seaseis -m rms_magg)
//
void params_mod_rms_magg_( csParamDef* pdef ) {
  pdef->setModule( "RMS_MAGG", "RMS amplitude in a window centred on the picked time",
                   "Stores the RMS of the samples between pick - length/2 and pick + length/2 in a trace header. Trace samples are not changed." );

  pdef->addDoc("Pick time [ms] is read from the trace header 'hdr_pick' (default time_pick, as written by module PICKING).");
  pdef->addDoc("Window samples: round((pick - length/2)/dt) ... round((pick + length/2)/dt), clipped to the trace.");
  pdef->addDoc("Traces without pick (pick <= 0, or window outside the trace) get 'null_value'. The log reports how many traces were affected.");

  pdef->addParam( "length", "Total window length [ms], centred on the pick", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Window length [ms]" );

  pdef->addParam( "hdr_pick", "Trace header containing the picked time [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "time_pick", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "hdr_rms", "Trace header where the RMS value is stored", NUM_VALUES_FIXED, "Created if it does not exist" );
  pdef->addValue( "rms_pick", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "remove_mean", "Remove the mean of the window before computing the RMS", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "RMS of the samples" );
  pdef->addOption( "yes", "RMS of the samples minus the window mean (standard deviation)" );

  pdef->addParam( "zero_pick", "How to treat pick time = 0", NUM_VALUES_FIXED );
  pdef->addValue( "ignore", VALTYPE_OPTION );
  pdef->addOption( "ignore", "Pick = 0 means 'no pick': RMS = null_value" );
  pdef->addOption( "use", "Pick = 0 is a valid time (window is clipped at the start of the trace)" );

  pdef->addParam( "null_value", "Value stored for traces without valid pick", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER );

  pdef->addParam( "hdr_nsamples", "Optional trace header for the number of samples used (less than full window if clipped)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );
}

//************************************************************************************************
// Start exec phase
//
bool start_exec_mod_rms_magg_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

//************************************************************************************************
// Cleanup phase
//
void cleanup_mod_rms_magg_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  writer->line( "RMS_MAGG: %d traces;  %d without pick;  %d with window outside trace;  %d with window clipped at trace start/end",
                vars->numTraces, vars->numNoPick, vars->numOutside, vars->numClipped );
  delete vars; vars = NULL;
}

extern "C" void _params_mod_rms_magg_( csParamDef* pdef ) {
  params_mod_rms_magg_( pdef );
}
extern "C" void _init_mod_rms_magg_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_rms_magg_( param, env, writer );
}
extern "C" bool _start_exec_mod_rms_magg_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_rms_magg_( env, writer );
}
extern "C" void _exec_mod_rms_magg_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_rms_magg_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_rms_magg_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_rms_magg_( env, writer );
}
