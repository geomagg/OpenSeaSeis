/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module SPK_DECON: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstring>
#include <cmath>
#include <vector>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_DECON
 *
 * Multi-window predictive deconvolution from SEISPAK process JHHDCON.
 * Parameters, defaults and checks follow the SEISPAK edit program (ejhhdcon.f);
 * computation is done by the Fortran routines in fortran/.
 */

extern "C" {
  void spkdcon_( float* s, int* nt, float* smp, int* nstm, float* stm, float* etm, float* vel,
                 float* flngth, float* gap, float* wnd, int* mode, float* fr, float* rng, int* irang,
                 float* r, float* y, float* x, float* f, float* aigap, float* aigsq, int* ncnt, int* ipflg );
}

namespace mod_spk_decon {
  static int const MAX_WIN = 10;
  struct VariableStruct {
    int nstm;
    int mode;
    float fr;
    float stm[MAX_WIN], etm[MAX_WIN], vel[MAX_WIN], flngth[MAX_WIN], gap[MAX_WIN];
    float wnd[2*MAX_WIN];
    float aigap[MAX_WIN], aigsq[MAX_WIN];
    int   ncnt[MAX_WIN];
    int   ipflg;
    int   hdrId_offset;
    std::vector<float> r, y, x, f;
  };
}
using mod_spk_decon::VariableStruct;
using mod_spk_decon::MAX_WIN;

static int readList( csParamManager* param, char const* name, float* values, csLogWriter* writer ) {
  if( !param->exists( name ) ) return 0;
  int n = param->getNumValues( name );
  if( n > 2*MAX_WIN ) writer->error("Too many values for parameter '%s' (%d)", name, n);
  for( int i = 0; i < n; i++ ) param->getFloat( name, &values[i], i );
  return n;
}

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_spk_decon_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  int   nt  = shdr->numSamples;
  float smp = shdr->sampleInt;
  float tmx = (nt-1)*smp;

  // Defaults as in ejhhdcon.f
  float stm[2*MAX_WIN], etm[2*MAX_WIN], vel[2*MAX_WIN], flngth[2*MAX_WIN], gap[2*MAX_WIN], wnd[4*MAX_WIN];
  for( int i = 0; i < 2*MAX_WIN; i++ ) {
    stm[i] = 0.0f; etm[i] = tmx; vel[i] = 0.0f; flngth[i] = -1.0f; gap[i] = 1.0f;
  }
  for( int i = 0; i < 4*MAX_WIN; i++ ) wnd[i] = -1.0f;
  float stblz = 0.1f;
  vars->mode  = 1;

  if( param->exists("mode") ) param->getInt( "mode", &vars->mode );
  if( vars->mode < 0 || vars->mode > 3 ) writer->error("Parameter 'mode' must be 0, 1, 2 or 3");
  if( param->exists("stabilize") ) param->getFloat( "stabilize", &stblz );

  int nstm = readList( param, "stm", stm, writer );
  if( nstm == 0 ) nstm = 1;
  int netm = readList( param, "etm", etm, writer );
  if( netm == 0 ) netm = 1;
  int nvel = readList( param, "vel", vel, writer );
  if( nvel == 0 ) nvel = 1;
  int nflngth = readList( param, "flngth", flngth, writer );
  int ngap = readList( param, "gap", gap, writer );
  if( ngap == 0 ) ngap = nstm;
  int nwnd = readList( param, "dwindow", wnd, writer );

  //------------------------------------------------------------
  // Checks as in ejhhdcon.f
  if( nstm > MAX_WIN ) writer->error("Maximum number of windows is %d", MAX_WIN);
  if( nstm != netm || nstm != nvel || nstm != nflngth || nstm != ngap || 2*nstm != nwnd ) {
    writer->error("Argument count error: DWINDOW must have twice as many values as STM, ETM, VEL, FLNGTH and GAP.\n"
                  "NSTM,NETM,NVEL: %d %d %d   NGAP,NWND,NFLNGTH: %d %d %d", nstm, netm, nvel, ngap, nwnd, nflngth);
  }
  if( etm[netm-1] > 0.0f && etm[netm-1] < tmx ) {
    writer->warning("Trace longer than last ETM (%f < %f ms)", etm[netm-1], tmx);
  }
  for( int i = 0; i < netm; i++ ) if( etm[i] <= 0.0f ) etm[i] = tmx;
  for( int i = 0; i < nstm; i++ ) {
    float w1 = wnd[2*i], w2 = wnd[2*i+1];
    if( w2 < w1 ) writer->error("Window %d: DWINDOW pairs must increase (%f %f)", i+1, w1, w2);
    if( w2 > tmx ) writer->error("Window %d: DWINDOW extends beyond trace (%f > %f ms)", i+1, w2, tmx);
    if( flngth[i] <= 0.0f ) writer->error("Window %d: FLNGTH must be > 0 ms", i+1);
    if( flngth[i] > (w2 - w1) ) writer->error("Window %d: filter must fit in design window (FLNGTH %f > %f ms)", i+1, flngth[i], w2-w1);
    if( etm[i] <= stm[i] ) writer->error("Window %d: STM > ETM", i+1);
    if( i > 0 && stm[i] > etm[i-1] ) writer->error("ETM must overlap subsequent STM: ETM(%d)=%f, STM(%d)=%f", i, etm[i-1], i+1, stm[i]);
    if( vars->mode == 2 && gap[i] >= flngth[i] ) writer->error("Window %d: GAP (%f ms) must be shorter than FLNGTH (%f ms)", i+1, gap[i], flngth[i]);
  }

  vars->nstm = nstm;
  float frx = 0.01f*stblz;
  vars->fr = 1.0f + frx*frx;         // as in JHHDCON: FR=0.01*STBLZ; FR=1+FR*FR
  int maxLngf = 1;
  for( int i = 0; i < nstm; i++ ) {
    vars->stm[i] = stm[i]; vars->etm[i] = etm[i]; vars->vel[i] = vel[i];
    vars->flngth[i] = flngth[i]; vars->gap[i] = gap[i];
    vars->wnd[2*i] = wnd[2*i]; vars->wnd[2*i+1] = wnd[2*i+1];
    vars->aigap[i] = 0.0f; vars->aigsq[i] = 0.0f; vars->ncnt[i] = 0;
    int lngf = (int)lrintf( flngth[i]/smp );
    if( lngf > maxLngf ) maxLngf = lngf;
  }
  if( maxLngf > 4000 ) writer->error("Filter too long (%d samples > 4000)", maxLngf);
  vars->ipflg = 1;

  //------------------------------------------------------------
  // Offset header (SEISPAK 'RANG'), used with VEL > 0
  vars->hdrId_offset = -1;
  std::string hdrName = "offset";
  if( param->exists("hdr_offset") ) param->getString( "hdr_offset", &hdrName );
  if( hdef->headerExists( hdrName ) ) {
    vars->hdrId_offset = hdef->headerIndex( hdrName );
  }
  else {
    bool useVel = false;
    for( int i = 0; i < nstm; i++ ) if( vel[i] > 0.0f ) useVel = true;
    if( useVel ) writer->warning("Offset header '%s' not found: windows will not be shifted with offset (VEL ignored)", hdrName.c_str());
  }

  vars->r.resize( nt + maxLngf + 1 );
  vars->y.resize( nt + 1 );
  vars->x.resize( (size_t)nt * nstm );
  vars->f.resize( maxLngf + 1 );

  //------------------------------------------------------------
  char const* modeText[4] = { "0 (output autocorrelation)", "1 (gap = GAP-th zero crossing)",
                              "2 (gap in ms)", "3 (gap = max between zero crossings)" };
  writer->line("  Mode:       %s", modeText[vars->mode]);
  writer->line("  Stabilize:  %f   (FR = %f)", stblz, vars->fr);
  writer->line("  Offset hdr: %s", vars->hdrId_offset >= 0 ? hdrName.c_str() : "(none)");
  writer->line("  Win     STM      ETM      VEL   FLNGTH      GAP   DWINDOW");
  for( int i = 0; i < nstm; i++ ) {
    writer->line("  %2d  %7.1f  %7.1f  %7.1f  %7.1f  %7.2f   %7.1f %7.1f", i+1,
                 stm[i], etm[i], vel[i], flngth[i], gap[i], wnd[2*i], wnd[2*i+1]);
  }
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_spk_decon_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;
  csTrace* trace = traceGather->trace(0);
  float* samples = trace->getTraceSamples();
  int   nt  = shdr->numSamples;
  float smp = shdr->sampleInt;

  float rng = 0.0f;
  int irang = 0;
  if( vars->hdrId_offset >= 0 ) {
    rng = (float)trace->getTraceHeader()->doubleValue( vars->hdrId_offset );
    irang = 1;
    if( rng != rng ) return;   // NaN offset: leave trace unchanged (as JHHDCON)
  }
  spkdcon_( samples, &nt, &smp, &vars->nstm, vars->stm, vars->etm, vars->vel, vars->flngth, vars->gap,
            vars->wnd, &vars->mode, &vars->fr, &rng, &irang,
            &vars->r[0], &vars->y[0], &vars->x[0], &vars->f[0],
            vars->aigap, vars->aigsq, vars->ncnt, &vars->ipflg );
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_spk_decon_( csParamDef* pdef ) {
  pdef->setModule( "SPK_DECON", "Multi-window predictive deconvolution (SEISPAK JHHDCON)",
                   "Wiener prediction-error filter designed in each design window (DWINDOW) and applied in the zone STM-ETM; "
                   "zones are merged with a linear taper where they overlap. With VEL > 0, windows are shifted by the hyperbolic "
                   "moveout of the offset header. Parameters, defaults and checks follow SEISPAK (ejhhdcon.f). "
                   "Up to 10 windows: STM, ETM, VEL, FLNGTH and GAP take one value per window, DWINDOW two values per window." );

  pdef->addParam( "mode", "Gap definition", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "0: output autocorrelation, 1: gap at GAP-th zero crossing of autocorrelation, "
                  "2: gap in ms, 3: gap at maximum of autocorrelation between zero crossings 2*GAP and 2*GAP+1" );

  pdef->addParam( "stabilize", "Stabilisation (white noise)", NUM_VALUES_FIXED, "FR = 1 + (0.01*STABILIZE)^2 multiplies the zero lag" );
  pdef->addValue( "0.1", VALTYPE_NUMBER, "Stabilisation" );

  pdef->addParam( "stm", "Start time of application zone, per window [ms]", NUM_VALUES_VARIABLE );
  pdef->addValue( "0", VALTYPE_NUMBER, "STM [ms]" );

  pdef->addParam( "etm", "End time of application zone, per window [ms]", NUM_VALUES_VARIABLE, "0 = end of trace" );
  pdef->addValue( "0", VALTYPE_NUMBER, "ETM [ms]" );

  pdef->addParam( "vel", "Velocity for moveout of windows, per window [m/s]", NUM_VALUES_VARIABLE, "0 = no moveout" );
  pdef->addValue( "0", VALTYPE_NUMBER, "VEL" );

  pdef->addParam( "flngth", "Filter length, per window [ms]", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "FLNGTH [ms]" );

  pdef->addParam( "gap", "Prediction gap, per window", NUM_VALUES_VARIABLE, "In ms for mode 2, number of zero crossings for modes 1 and 3" );
  pdef->addValue( "1", VALTYPE_NUMBER, "GAP" );

  pdef->addParam( "dwindow", "Design window, two values per window [ms]", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "Start [ms]" );
  pdef->addValue( "", VALTYPE_NUMBER, "End [ms]" );

  pdef->addParam( "hdr_offset", "Trace header with offset (SEISPAK RANG)", NUM_VALUES_FIXED, "Used only with VEL > 0" );
  pdef->addValue( "offset", VALTYPE_STRING, "Header name" );
}

//************************************************************************************************
// Start exec phase
//*************************************************************************************************
bool start_exec_mod_spk_decon_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

//************************************************************************************************
// Cleanup phase: gap statistics (as printed by JHHDCON after the last trace)
//*************************************************************************************************
void cleanup_mod_spk_decon_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  float smp = env->superHeader->sampleInt;
  if( vars->nstm > 0 && vars->mode != 0 ) {
    writer->line("");
    writer->line("  Window   Average gap [ms]   Std. dev. [ms]   Traces");
    for( int i = 0; i < vars->nstm; i++ ) {
      if( vars->ncnt[i] <= 0 ) continue;
      double avg = vars->aigap[i] / vars->ncnt[i];
      double sq  = vars->aigsq[i] / vars->ncnt[i] - avg*avg;
      double sdv = ( sq > 0.0 ) ? sqrt( sq ) : 0.0;
      writer->line("  %4d     %10.1f         %10.1f     %8d", i+1, avg*smp, sdv*smp, vars->ncnt[i]);
    }
  }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_spk_decon_( csParamDef* pdef ) {
  params_mod_spk_decon_( pdef );
}
extern "C" void _init_mod_spk_decon_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spk_decon_( param, env, writer );
}
extern "C" bool _start_exec_mod_spk_decon_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spk_decon_( env, writer );
}
extern "C" void _exec_mod_spk_decon_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_decon_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_decon_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spk_decon_( env, writer );
}
