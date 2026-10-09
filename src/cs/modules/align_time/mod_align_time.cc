/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/*
 * Module ALIGN_TIME: aligns traces in absolute time (headers time_samp1 + time_samp1_us, or
 * time_year/day/hour/min/sec/msec/usec).
 * Author: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026.
 * Same idea as the SeaView plugin "Alinhar traços no tempo absoluto" (csProcessingAlignTime).
 *
 * mode ensemble: inside each ensemble (or the whole input, without ensembles) every trace is shifted so
 *                that sample 0 is the same instant: the LATEST start time of the ensemble (or the earliest).
 *                Optionally samples after the earliest end are zeroed (cut yes).
 * mode grid:     each trace is shifted on its own so that it starts at the next multiple of 'grid' [ms]
 *                counted from 1970-01-01 00:00:00 UTC. Traces from different nodes then share exactly the
 *                same sample times, trace by trace, without buffering (good for long records).
 * Shift: new(k) = old(k + d), d = (t_new - t_old)/dt samples (linear interpolation, or nearest sample).
 * time_samp1/time_samp1_us are set to the new start time. Samples without input are zero.
 */

#include "cseis_includes.h"
#include "csStandardHeaders.h"
#include <cmath>
#include <cstring>
#include <ctime>
#include <vector>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

namespace mod_align_time {
  struct VariableStruct {
    int mode;                 // 1 ensemble, 2 grid
    bool latest;              // ensemble: common start = latest (true) or earliest start
    bool cut;                 // ensemble: zero after the earliest end
    bool linear;              // linear interpolation (else nearest sample)
    double grid;              // [s]
    int hdrId_t1, hdrId_t1us;
    int hdrId_year, hdrId_day, hdrId_hour, hdrId_min, hdrId_sec, hdrId_msec, hdrId_usec;
    int hdrId_shift;          // header 'align_shift_ms' (applied shift), -1 if not requested
    long numTraces, numNoTime, numShifted;
    double maxShift;          // [s]
  };
  static int const MODE_ENSEMBLE = 1;
  static int const MODE_GRID     = 2;
  static double startTime( csTraceHeader const* h, VariableStruct const* v );
  static void setStartTime( csTraceHeader* h, VariableStruct const* v, double t );
  static void shiftTrace( float* s, int ns, double d, bool linear, std::vector<float>& work );
}
using namespace mod_align_time;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_align_time_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csTraceHeaderDef* hdef = env->headerDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars   = new VariableStruct();
  edef->setVariables( vars );

  vars->mode = MODE_ENSEMBLE; vars->latest = true; vars->cut = true; vars->linear = true;
  vars->grid = shdr->sampleInt / 1000.0;
  vars->numTraces = vars->numNoTime = vars->numShifted = 0;
  vars->maxShift = 0.0;
  string text;

  if( param->exists( "mode" ) ) {
    param->getString( "mode", &text );
    if( !text.compare( "ensemble" ) ) vars->mode = MODE_ENSEMBLE;
    else if( !text.compare( "grid" ) ) vars->mode = MODE_GRID;
    else writer->error( "Unknown option for 'mode': '%s'", text.c_str() );
  }
  if( param->exists( "start" ) ) {
    param->getString( "start", &text );
    if( !text.compare( "latest" ) ) vars->latest = true;
    else if( !text.compare( "earliest" ) ) vars->latest = false;
    else writer->error( "Unknown option for 'start': '%s'", text.c_str() );
  }
  if( param->exists( "cut" ) ) {
    param->getString( "cut", &text );
    vars->cut = !text.compare( "yes" );
  }
  if( param->exists( "interp" ) ) {
    param->getString( "interp", &text );
    if( !text.compare( "linear" ) ) vars->linear = true;
    else if( !text.compare( "nearest" ) ) vars->linear = false;
    else writer->error( "Unknown option for 'interp': '%s'", text.c_str() );
  }
  if( param->exists( "grid" ) ) {
    double g;
    param->getDouble( "grid", &g );
    if( g <= 0 ) writer->error( "grid must be > 0" );
    vars->grid = g / 1000.0;
  }

  if( vars->mode == MODE_ENSEMBLE ) edef->setTraceSelectionMode( TRCMODE_ENSEMBLE );
  else edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->hdrId_t1   = hdef->headerExists( HDR_TIME_SAMP1.name ) ? hdef->headerIndex( HDR_TIME_SAMP1.name ) : -1;
  vars->hdrId_t1us = hdef->headerExists( HDR_TIME_SAMP1_US.name ) ? hdef->headerIndex( HDR_TIME_SAMP1_US.name ) : -1;
  vars->hdrId_year = hdef->headerExists( "time_year" ) ? hdef->headerIndex( "time_year" ) : -1;
  vars->hdrId_day  = hdef->headerExists( "time_day" )  ? hdef->headerIndex( "time_day" )  : -1;
  vars->hdrId_hour = hdef->headerExists( "time_hour" ) ? hdef->headerIndex( "time_hour" ) : -1;
  vars->hdrId_min  = hdef->headerExists( "time_min" )  ? hdef->headerIndex( "time_min" )  : -1;
  vars->hdrId_sec  = hdef->headerExists( "time_sec" )  ? hdef->headerIndex( "time_sec" )  : -1;
  vars->hdrId_msec = hdef->headerExists( "time_msec" ) ? hdef->headerIndex( "time_msec" ) : -1;
  vars->hdrId_usec = hdef->headerExists( "time_usec" ) ? hdef->headerIndex( "time_usec" ) : -1;
  if( vars->hdrId_t1 < 0 && vars->hdrId_day < 0 && vars->hdrId_hour < 0 && vars->hdrId_min < 0 && vars->hdrId_sec < 0 ) {
    writer->error( "No time headers: time_samp1 (+time_samp1_us) or time_year/day/hour/min/sec are required" );
  }
  // The new start time is always written to time_samp1 / time_samp1_us
  if( vars->hdrId_t1 < 0 )   vars->hdrId_t1   = hdef->addStandardHeader( HDR_TIME_SAMP1.name );
  if( vars->hdrId_t1us < 0 ) vars->hdrId_t1us = hdef->addStandardHeader( HDR_TIME_SAMP1_US.name );
  if( !hdef->headerExists( "align_shift_ms" ) ) hdef->addHeader( TYPE_FLOAT, "align_shift_ms", "Time shift applied by ALIGN_TIME [ms]: new start - old start" );
  vars->hdrId_shift = hdef->headerIndex( "align_shift_ms" );

  if( vars->mode == MODE_ENSEMBLE ) {
    writer->line( "  Mode ensemble: common start = %s start of each ensemble; cut at the earliest end: %s; interpolation: %s",
                  vars->latest ? "latest" : "earliest", vars->cut ? "yes" : "no", vars->linear ? "linear" : "nearest sample" );
  }
  else {
    writer->line( "  Mode grid: each trace starts at the next multiple of %g ms since 1970-01-01 UTC; interpolation: %s",
                  vars->grid * 1000.0, vars->linear ? "linear" : "nearest sample" );
  }
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_align_time_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader const* shdr = env->superHeader;
  int ntr = traceGather->numTraces();
  int ns  = shdr->numSamples;
  double dt = shdr->sampleInt / 1000.0;
  std::vector<float> work( ns );

  std::vector<double> t0( ntr );
  for( int i = 0; i < ntr; i++ ) t0[i] = startTime( traceGather->trace(i)->getTraceHeader(), vars );

  double tNew = 0.0, tEnd = 0.0;
  if( vars->mode == MODE_ENSEMBLE ) {
    bool any = false;
    for( int i = 0; i < ntr; i++ ) {
      if( std::isnan( t0[i] ) ) continue;
      double te = t0[i] + ( ns - 1 ) * dt;
      if( !any ) { tNew = t0[i]; tEnd = te; any = true; }
      else {
        tNew = vars->latest ? std::max( tNew, t0[i] ) : std::min( tNew, t0[i] );
        tEnd = std::min( tEnd, te );
      }
    }
    if( !any ) { vars->numNoTime += ntr; vars->numTraces += ntr; return; }
    if( tEnd <= tNew ) {
      writer->error( "Ensemble with %d traces: the traces do not overlap in time (start %.3f s after the earliest end).\n"
                     "The ensemble probably contains several files/time chunks: define one ensemble per file, e.g. ENS_DEFINE header fileno",
                     ntr, tNew - tEnd );
    }
  }

  for( int i = 0; i < ntr; i++ ) {
    vars->numTraces++;
    if( std::isnan( t0[i] ) ) { vars->numNoTime++; continue; }
    double tn = tNew;
    if( vars->mode == MODE_GRID ) {
      // next multiple of the grid (a start already on the grid, within 1 us, is kept)
      double k = ceil( t0[i] / vars->grid - 1.0e-6 / vars->grid );
      tn = k * vars->grid;
    }
    double d = ( tn - t0[i] ) / dt;            // samples
    float* s = traceGather->trace(i)->getTraceSamples();
    if( fabs( d ) > 1.0e-6 ) {
      shiftTrace( s, ns, d, vars->linear, work );
      vars->numShifted++;
    }
    if( vars->mode == MODE_ENSEMBLE && vars->cut ) {
      int nKeep = ( tEnd > tn ) ? (int)floor( ( tEnd - tn ) / dt + 1.0e-6 ) + 1 : 0;
      for( int k = std::max( 0, nKeep ); k < ns; k++ ) s[k] = 0.0f;
    }
    csTraceHeader* h = traceGather->trace(i)->getTraceHeader();
    setStartTime( h, vars, tn );
    h->setFloatValue( vars->hdrId_shift, (float)( ( tn - t0[i] ) * 1000.0 ) );
    vars->maxShift = std::max( vars->maxShift, fabs( tn - t0[i] ) );
  }
  (void)edef;
}

namespace mod_align_time {
/** Absolute start time [s since 1970], NaN if the trace has no time */
static double startTime( csTraceHeader const* h, VariableStruct const* v ) {
  if( v->hdrId_t1 >= 0 && h->intValue( v->hdrId_t1 ) != 0 ) {
    return (double)h->intValue( v->hdrId_t1 ) + ( v->hdrId_t1us >= 0 ? 1.0e-6 * h->intValue( v->hdrId_t1us ) : 0.0 );
  }
  if( v->hdrId_day < 0 && v->hdrId_hour < 0 && v->hdrId_min < 0 && v->hdrId_sec < 0 ) return NAN;
  int year = ( v->hdrId_year >= 0 ) ? h->intValue( v->hdrId_year ) : 1970;
  if( year <= 0 ) year = 1970;
  else if( year < 100 ) year += 2000;
  int day = ( v->hdrId_day >= 0 ) ? std::max( 1, h->intValue( v->hdrId_day ) ) : 1;
  // days from 1970-01-01 to year-01-01
  long days = 0;
  for( int y = 1970; y < year; y++ ) days += ( ( y % 4 == 0 && y % 100 != 0 ) || y % 400 == 0 ) ? 366 : 365;
  double t = 86400.0 * ( days + day - 1 );
  if( v->hdrId_hour >= 0 ) t += 3600.0 * h->doubleValue( v->hdrId_hour );
  if( v->hdrId_min >= 0 )  t += 60.0 * h->doubleValue( v->hdrId_min );
  if( v->hdrId_sec >= 0 )  t += h->doubleValue( v->hdrId_sec );
  if( v->hdrId_msec >= 0 ) t += 1.0e-3 * h->doubleValue( v->hdrId_msec );
  if( v->hdrId_usec >= 0 ) t += 1.0e-6 * h->doubleValue( v->hdrId_usec );
  return t;
}
static void setStartTime( csTraceHeader* h, VariableStruct const* v, double t ) {
  double s = floor( t );
  int us = (int)lround( ( t - s ) * 1.0e6 );
  if( us >= 1000000 ) { s += 1.0; us -= 1000000; }
  h->setIntValue( v->hdrId_t1, (int)s );
  h->setIntValue( v->hdrId_t1us, us );
}
/** new(k) = old(k + d) */
static void shiftTrace( float* s, int ns, double d, bool linear, std::vector<float>& work ) {
  if( !linear ) d = (double)lround( d );
  int id = (int)floor( d );
  double w = d - id;
  if( w < 1.0e-6 ) w = 0.0;
  for( int k = 0; k < ns; k++ ) {
    int j = k + id;                          // samples outside the trace are zero
    double a0 = ( j >= 0 && j < ns ) ? s[j] : 0.0;
    double a1 = ( j + 1 >= 0 && j + 1 < ns ) ? s[j+1] : 0.0;
    work[k] = (float)( ( w > 0.0 ) ? ( 1.0 - w ) * a0 + w * a1 : a0 );
  }
  memcpy( s, &work[0], ns * sizeof(float) );
}
}

//*************************************************************************************************
// Parameter definition
//*************************************************************************************************
void params_mod_align_time_( csParamDef* pdef ) {
  pdef->setModule( "ALIGN_TIME", "Align traces in absolute time",
    "Shifts the samples so that traces share the same absolute sample times (headers time_samp1/time_samp1_us, "
    "or time_year/day/hour/min/sec). The new start time is written to time_samp1/time_samp1_us and the applied shift "
    "to header align_shift_ms. Same idea as the SeaView plugin 'Alinhar traços no tempo absoluto'." );

  pdef->addParam( "mode", "Alignment mode", NUM_VALUES_FIXED );
  pdef->addValue( "ensemble", VALTYPE_OPTION );
  pdef->addOption( "ensemble", "All traces of each ensemble (or of the whole input, without ensembles) start at the same instant" );
  pdef->addOption( "grid", "Each trace starts at the next multiple of 'grid' since 1970-01-01 UTC (trace by trace, no buffering)" );

  pdef->addParam( "start", "Common start time (mode ensemble)", NUM_VALUES_FIXED );
  pdef->addValue( "latest", VALTYPE_OPTION );
  pdef->addOption( "latest", "Latest start of the ensemble: no trace needs samples before its first sample" );
  pdef->addOption( "earliest", "Earliest start of the ensemble: later traces are padded with zeros at the start" );

  pdef->addParam( "cut", "Zero samples after the earliest end (mode ensemble)", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "All traces end where the earliest-ending trace ends (samples after it are zero)" );
  pdef->addOption( "no", "Keep the samples after the common end" );

  pdef->addParam( "grid", "Time grid [ms] (mode grid)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Grid step [ms]. Default: sample interval. Example: 10000 for 10-s files" );

  pdef->addParam( "interp", "Interpolation of fractional shifts", NUM_VALUES_FIXED );
  pdef->addValue( "linear", VALTYPE_OPTION );
  pdef->addOption( "linear", "Linear interpolation between samples" );
  pdef->addOption( "nearest", "Shift by the nearest whole number of samples (start time is still set exactly)" );
}

bool start_exec_mod_align_time_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

void cleanup_mod_align_time_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  writer->line( "ALIGN_TIME: %ld traces, %ld shifted, largest shift %.3f ms", vars->numTraces, vars->numShifted, vars->maxShift * 1000.0 );
  if( vars->numNoTime > 0 ) writer->line( "  %ld trace(s) without time headers: not changed", vars->numNoTime );
  delete vars;
}

extern "C" void _params_mod_align_time_( csParamDef* pdef ) {
  params_mod_align_time_( pdef );
}
extern "C" void _init_mod_align_time_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_align_time_( param, env, writer );
}
extern "C" bool _start_exec_mod_align_time_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_align_time_( env, writer );
}
extern "C" void _exec_mod_align_time_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_align_time_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_align_time_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_align_time_( env, writer );
}
