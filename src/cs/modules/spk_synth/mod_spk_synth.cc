/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module SPK_SYNTH: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_SYNTH
 *
 * Input module: synthetic data of SEISPAK process SYNTH. Points, lines, hyperbolas,
 * hyperboloids (HYPERBOLA3D) and DMO bin responses (BIN3D) are placed on HMR records of HMT
 * traces, as spikes or gaussians. With VFILE an ASCII velocity file is converted to traces.
 * Computation is done by the original Fortran routines; the output data set (file ONLINE1 in
 * SEISPAK) is kept in memory and accessed through GETTRC/PUTTRC.
 */

extern "C" {
  void spksynth_( float* s, char* vfile, float* gssn, float* ggd, float* spj, float* spi, float* zrngtrc, float* tpwr,
                  int* lrm, float* smp, int* ihmr, int* ihmt,
                  int* mpoint, float* point, int* mline, float* rline, int* mhypb, float* hypb,
                  int* mhypb3d, float* hypb3d, int* mbin3d, float* bin3d,
                  int* nrec, int* ierr, size_t vfileLen );
}

namespace mod_spk_synth {
  static int const LHDR = 7;            // RNUM TNUM RANG SX GX ILIN XLIN
  static char const* SPK_NAMES[LHDR] = { "rnum", "tnum", "rang", "sx", "gx", "ilin", "xlin" };

  // Data set of SYNTH (ONLINE1): records x traces x (LHDR + NT), grown on demand (VFILE mode)
  struct DataSet {
    int hmt;
    int ns;                              // LHDR + NT
    std::vector< std::vector<float> > rec;
    int numOutside;
  };
  static DataSet* g_data = NULL;

  struct VariableStruct {
    // Parameters (SEISPAK names)
    int   hmr, hmt, lrm, nt;
    float smp, ggd, gssn, tpwr, spj, spi, zrngtrc;
    std::string vfile;
    std::vector<float> point, line, hyperbola, hyperbola3d, bin3d;
    // State
    DataSet data;
    int nrec;
    bool computed, atEOF;
    int traceCounter;
    // Headers
    int hdrId_trcno, hdrId_ffid, hdrId_chan, hdrId_offset, hdrId_soux, hdrId_recx, hdrId_binx, hdrId_biny, hdrId_row, hdrId_col;
    int hdrId_spk[LHDR];
  };

  // Reads all lines of a parameter: each line is one group (one SEISPAK keyword) with a multiple of
  // 'nset' values (after 'nfirst' leading values). Appends N, v1..vN per line as ESYNTH does.
  static void readGroups( csParamManager* param, csLogWriter* writer, char const* name, int nfirst, int nset, std::vector<float>& out ) {
    out.clear();
    int nLines = param->getNumLines( name );
    for( int iline = 0; iline < nLines; iline++ ) {
      int n = param->getNumValues( name, iline );
      if( n < nfirst + nset || ( n - nfirst ) % nset != 0 ) {
        if( nfirst > 0 ) writer->error("Parameter '%s', line %d: %d values given. Expected %d + a multiple of %d.", name, iline+1, n, nfirst, nset);
        else writer->error("Parameter '%s', line %d: %d values given. Expected a multiple of %d.", name, iline+1, n, nset);
      }
      out.push_back( (float)n );
      for( int i = 0; i < n; i++ ) {
        float v;
        param->getFloatAtLine( name, &v, iline, i );
        out.push_back( v );
      }
    }
  }
}
using namespace mod_spk_synth;

//--------------------------------------------------------------------------------
// Called by the Fortran routines (original SEISPAK I/O on the ONLINE1 data set)
//
extern "C" void gettrc_( int* ir, int* it, float* s ) {
  DataSet* d = g_data;
  int irec = *ir, itrc = *it;
  if( irec < 1 || itrc < 1 || itrc > d->hmt ) {
    d->numOutside++;
    memset( s, 0, d->ns * sizeof(float) );
    return;
  }
  if( irec > (int)d->rec.size() ) {
    memset( s, 0, d->ns * sizeof(float) );
    s[0] = (float)irec;
    s[1] = (float)itrc;
    return;
  }
  memcpy( s, &d->rec[irec-1][ (size_t)(itrc-1) * d->ns ], d->ns * sizeof(float) );
}

extern "C" void puttrc_( int* ir, int* it, float* s ) {
  DataSet* d = g_data;
  int irec = *ir, itrc = *it;
  if( irec < 1 || itrc < 1 || itrc > d->hmt ) {
    d->numOutside++;
    return;
  }
  while( (int)d->rec.size() < irec ) d->rec.push_back( std::vector<float>( (size_t)d->hmt * d->ns, 0.0f ) );
  memcpy( &d->rec[irec-1][ (size_t)(itrc-1) * d->ns ], s, d->ns * sizeof(float) );
}

//*************************************************************************************************
// Init phase
//
void init_mod_spk_synth_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csTraceHeaderDef* hdef = env->headerDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  edef->setExecType( EXEC_TYPE_INPUT );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->computed = false;
  vars->atEOF    = false;
  vars->traceCounter = 0;
  vars->nrec = 0;

  // Presets as in ESYNTH
  vars->hmr = 0; vars->hmt = 0; vars->lrm = 0;
  vars->smp = 0; vars->ggd = 0;
  vars->gssn = 0; vars->tpwr = 0; vars->spj = 1; vars->spi = 0; vars->zrngtrc = 1;
  vars->vfile = "NONE";

  if( param->exists("vfile") ) {
    param->getString( "vfile", &vars->vfile );
    if( vars->vfile.size() > 240 ) writer->error("VFILE name longer than 240 characters");
  }
  bool isVfile = vars->vfile.compare("NONE") != 0;

  float val;
  if( param->exists("hmr") ) { param->getFloat( "hmr", &val ); vars->hmr = (int)val; }
  if( param->exists("hmt") ) { param->getFloat( "hmt", &val ); vars->hmt = (int)val; }
  if( param->exists("lrm") ) { param->getFloat( "lrm", &val ); vars->lrm = (int)val; }
  if( param->exists("smp") ) param->getFloat( "smp", &vars->smp );
  if( param->exists("ggd") ) param->getFloat( "ggd", &vars->ggd );
  if( param->exists("gaussian") ) param->getFloat( "gaussian", &vars->gssn );
  if( param->exists("tpower") )   param->getFloat( "tpower", &vars->tpwr );
  if( param->exists("spj") )      param->getFloat( "spj", &vars->spj );
  if( param->exists("spi") )      param->getFloat( "spi", &vars->spi );
  if( param->exists("zrngtrc") )  param->getFloat( "zrngtrc", &vars->zrngtrc );

  readGroups( param, writer, "point",       0, 4, vars->point );
  readGroups( param, writer, "line",        0, 7, vars->line );
  readGroups( param, writer, "hyperbola",   0, 5, vars->hyperbola );
  readGroups( param, writer, "hyperbola3d", 0, 5, vars->hyperbola3d );
  readGroups( param, writer, "bin3d",       1, 7, vars->bin3d );

  // Checks as in ESYNTH
  if( !isVfile ) {
    if( vars->hmt <= 0 ) writer->error("HMT must be given > 0");
    if( vars->hmr <= 0 ) writer->error("HMR must be given > 0");
    if( vars->ggd <= 0 ) writer->error("GGD must be given > 0");
  }
  else {
    vars->hmt = 1;
    vars->hmr = 10000;
    vars->ggd = 1.0e30f;
    if( !vars->point.empty() || !vars->line.empty() || !vars->hyperbola.empty() || !vars->hyperbola3d.empty() || !vars->bin3d.empty() ) {
      writer->warning("VFILE given: POINT, LINE, HYPERBOLA, HYPERBOLA3D and BIN3D are ignored.");
    }
  }
  if( vars->smp <= 0 ) writer->error("SMP must be given > 0");
  if( vars->lrm <= 0 ) writer->error("LRM must be given > 0");
  vars->nt = 1 + (int)( vars->lrm / vars->smp );   // NT=1+(LRM/SMP) in SYNTH

  vars->data.hmt = vars->hmt;
  vars->data.ns  = LHDR + vars->nt;
  vars->data.numOutside = 0;

  shdr->numSamples = vars->nt;
  shdr->sampleInt  = vars->smp;
  shdr->domain     = DOMAIN_XT;

  if( isVfile ) {
    writer->line("  VFILE: %s  (one trace per record, %d samples, SMP %g ms)", vars->vfile.c_str(), vars->nt, vars->smp);
  }
  else {
    writer->line("  HMR %d  HMT %d  LRM %d ms  SMP %g ms  (%d samples)  GGD %g  GAUSSIAN %g  TPOWER %g  SPJ %g  SPI %g  ZRNGTRC %g",
                 vars->hmr, vars->hmt, vars->lrm, vars->smp, vars->nt, vars->ggd, vars->gssn, vars->tpwr, vars->spj, vars->spi, vars->zrngtrc);
    writer->line("  Values: POINT %d  LINE %d  HYPERBOLA %d  HYPERBOLA3D %d  BIN3D %d",
                 (int)vars->point.size(), (int)vars->line.size(), (int)vars->hyperbola.size(), (int)vars->hyperbola3d.size(), (int)vars->bin3d.size());
  }

  // Headers: SEISPAK global headers of ONLINE1 (spk_*, re-written by OUTPUT_SEISPAK) + SeaSeis standard headers
  for( int i = 0; i < LHDR; i++ ) {
    std::string name = std::string("spk_") + SPK_NAMES[i];
    if( !hdef->headerExists( name ) ) hdef->addHeader( TYPE_FLOAT, name, std::string("SEISPAK header ") + SPK_NAMES[i] + " (SYNTH)" );
    vars->hdrId_spk[i] = hdef->headerIndex( name );
  }
  vars->hdrId_trcno  = hdef->addStandardHeader( HDR_TRCNO.name );
  vars->hdrId_ffid   = hdef->addStandardHeader( HDR_FFID.name );
  vars->hdrId_chan   = hdef->addStandardHeader( HDR_CHAN.name );
  vars->hdrId_row    = hdef->addStandardHeader( HDR_ROW.name );
  vars->hdrId_col    = hdef->addStandardHeader( HDR_COL.name );
  vars->hdrId_offset = vars->hdrId_soux = vars->hdrId_recx = vars->hdrId_binx = vars->hdrId_biny = -1;
  if( isVfile ) {
    vars->hdrId_binx = hdef->addStandardHeader( HDR_BIN_X.name );
    vars->hdrId_biny = hdef->addStandardHeader( HDR_BIN_Y.name );
  }
  else {
    vars->hdrId_offset = hdef->addStandardHeader( HDR_OFFSET.name );
    vars->hdrId_soux   = hdef->addStandardHeader( HDR_SOU_X.name );
    vars->hdrId_recx   = hdef->addStandardHeader( HDR_REC_X.name );
  }
}

//*************************************************************************************************
// Exec phase
//
void exec_mod_spk_synth_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( vars->atEOF ) {
    traceGather->freeAllTraces();
    return;
  }

  if( !vars->computed ) {
    // Arrays as built by ESYNTH, plus the final zero added by SYNTH (POINT(MPOINT+1)=0 ...)
    std::vector<float> point = vars->point, line = vars->line, hyp = vars->hyperbola, hyp3d = vars->hyperbola3d, bin3d = vars->bin3d;
    int mpoint = (int)point.size(), mline = (int)line.size(), mhyp = (int)hyp.size(), mhyp3d = (int)hyp3d.size(), mbin3d = (int)bin3d.size();
    point.push_back(0); line.push_back(0); hyp.push_back(0); hyp3d.push_back(0); bin3d.push_back(0);

    char vfile[241];
    memset( vfile, ' ', 240 );
    vfile[240] = '\0';
    memcpy( vfile, vars->vfile.c_str(), vars->vfile.size() );

    std::vector<float> s( vars->data.ns, 0.0f );
    int nrec = 0, ierr = 0;
    g_data = &vars->data;
    spksynth_( &s[0], vfile, &vars->gssn, &vars->ggd, &vars->spj, &vars->spi, &vars->zrngtrc, &vars->tpwr,
               &vars->lrm, &vars->smp, &vars->hmr, &vars->hmt,
               &mpoint, &point[0], &mline, &line[0], &mhyp, &hyp[0], &mhyp3d, &hyp3d[0], &mbin3d, &bin3d[0],
               &nrec, &ierr, (size_t)240 );
    g_data = NULL;
    if( ierr == 1 ) writer->error("BIN3D: more records (HMR) than ranges given");
    if( ierr == 2 ) writer->error("BIN3D: range not found among the records (see standard output)");
    if( vars->vfile.compare("NONE") != 0 ) {
      if( nrec <= 0 ) writer->error("No velocity function read from VFILE '%s' (see standard output)", vars->vfile.c_str());
      vars->nrec = nrec;
    }
    else {
      vars->nrec = vars->hmr;
    }
    // Records never written (possible in VFILE mode) are empty
    while( (int)vars->data.rec.size() < vars->nrec ) vars->data.rec.push_back( std::vector<float>( (size_t)vars->hmt * vars->data.ns, 0.0f ) );
    if( vars->data.numOutside > 0 ) {
      writer->line("  NOTE: %d accesses to record/trace numbers outside the data set were ignored", vars->data.numOutside);
    }
    writer->line("  Generated %d records x %d traces", vars->nrec, vars->hmt);
    vars->computed = true;
  }

  int k = vars->traceCounter;
  int ntr = vars->nrec * vars->hmt;
  if( k >= ntr ) {
    vars->atEOF = true;
    traceGather->freeAllTraces();
    return;
  }
  if( k == ntr - 1 ) vars->atEOF = true;
  int irec = k / vars->hmt;
  int itrc = k % vars->hmt;
  float const* s = &vars->data.rec[irec][ (size_t)itrc * vars->data.ns ];

  csTrace* trace = traceGather->trace(0);
  memcpy( trace->getTraceSamples(), s + LHDR, vars->nt * sizeof(float) );
  csTraceHeader* trcHdr = trace->getTraceHeader();
  for( int i = 0; i < LHDR; i++ ) trcHdr->setFloatValue( vars->hdrId_spk[i], s[i] );
  trcHdr->setIntValue( vars->hdrId_trcno, k + 1 );
  trcHdr->setIntValue( vars->hdrId_ffid, (int)lround( s[0] ) );
  trcHdr->setIntValue( vars->hdrId_chan, (int)lround( s[1] ) );
  trcHdr->setIntValue( vars->hdrId_row,  (int)lround( s[5] ) );
  trcHdr->setIntValue( vars->hdrId_col,  (int)lround( s[6] ) );
  if( vars->hdrId_offset >= 0 ) {
    trcHdr->setFloatValue( vars->hdrId_offset, s[2] );
    trcHdr->setDoubleValue( vars->hdrId_soux, s[3] );
    trcHdr->setDoubleValue( vars->hdrId_recx, s[4] );
  }
  else {
    trcHdr->setDoubleValue( vars->hdrId_binx, s[3] );
    trcHdr->setDoubleValue( vars->hdrId_biny, s[4] );
  }
  vars->traceCounter += 1;
  if( vars->atEOF ) {
    vars->data.rec.clear();
  }
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_spk_synth_( csParamDef* pdef ) {
  pdef->setModule( "SPK_SYNTH", "Synthetic data: points, lines, hyperbolas (SEISPAK SYNTH)",
                   "Input module. Generates HMR records of HMT traces with points, lines, hyperbolas (one-way 'direct arrival' "
                   "hyperbolas), hyperboloids and DMO bin responses, as spikes or gaussians. With VFILE an ASCII velocity file "
                   "is converted to traces. Computation is done by the original Fortran routines." );
  pdef->addDoc("Each line of POINT, LINE, HYPERBOLA, HYPERBOLA3D or BIN3D is one SEISPAK keyword: several value sets on the same line, "
               "starting on different records, are interpolated over the records between them (SEISPAK continuation lines).");
  pdef->addDoc("Example (SEISPAK 'line 1 1 200 100 25 800 100' + continuation '3 3 300 50 20 600 150'):  line  1 1 200 100 25 800 100  3 3 300 50 20 600 150");
  pdef->addDoc("Headers: spk_rnum, spk_tnum, spk_rang, spk_sx, spk_gx, spk_ilin, spk_xlin (as ONLINE1; written back by OUTPUT_SEISPAK), "
               "and ffid = RNUM, chan = TNUM, offset = RANG, sou_x = SX, rec_x = GX, row = ILIN, col = XLIN (VFILE: bin_x, bin_y).");
  pdef->addDoc("RANG = (trace - ZRNGTRC)*GGD,  SX = SPI + (record-1)*SPJ*GGD,  GX = SX + RANG.");

  pdef->addParam( "hmr", "How many records to generate", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "HMR" );
  pdef->addParam( "hmt", "How many traces per record to generate", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "HMT" );
  pdef->addParam( "lrm", "Last millisecond on record", NUM_VALUES_FIXED, "Number of samples = 1 + LRM/SMP" );
  pdef->addValue( "", VALTYPE_NUMBER, "LRM [ms]" );
  pdef->addParam( "smp", "Sample interval in milliseconds", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "SMP [ms]" );
  pdef->addParam( "ggd", "Group spacing within each record", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "GGD" );

  pdef->addParam( "point", "Point(s): Rec Trc Time Amp", NUM_VALUES_VARIABLE, "4 values per set; several sets on one line are interpolated over the records between" );
  pdef->addValue( "", VALTYPE_NUMBER, "Rec, Trc, Time [ms], Amp, (Rec, Trc, Time, Amp, ...)" );
  pdef->addParam( "line", "Line(s): Rec Trc1 Time1 Amp1 Trc2 Time2 Amp2", NUM_VALUES_VARIABLE, "7 values per set; several sets on one line are interpolated over the records between" );
  pdef->addValue( "", VALTYPE_NUMBER, "Rec, Trc1, Time1 [ms], Amp1, Trc2, Time2 [ms], Amp2, ..." );
  pdef->addParam( "hyperbola", "Hyperbola(s): Rec Trc Time Amp Vel", NUM_VALUES_VARIABLE,
                  "One-way time t = sqrt(Time^2 + x^2/Vel^2), x = GGD*(trace - Trc). For a diffraction on a stacked section halve the velocity. "
                  "5 values per set; several sets on one line are interpolated over the records between" );
  pdef->addValue( "", VALTYPE_NUMBER, "Rec, Trc, Time [ms], Amp, Vel, ..." );
  pdef->addParam( "hyperbola3d", "Hyperboloid(s): Rec Trc Time Amp Vel", NUM_VALUES_VARIABLE,
                  "Second spatial dimension along the records, with record spacing = GGD. 5 values per set" );
  pdef->addValue( "", VALTYPE_NUMBER, "Rec, Trc, Time [ms], Amp, Vel, ..." );
  pdef->addParam( "bin3d", "DMO bin response: Range, then sets of trc vel z th ph sd sa", NUM_VALUES_VARIABLE,
                  "One line per common-range record (range = first value), followed by 7 values per bin: trace of the output bin, velocity, "
                  "perpendicular distance bin-reflector, angle z axis-reflector normal [deg], azimuth of the normal [deg], shot distance and azimuth from the bin" );
  pdef->addValue( "", VALTYPE_NUMBER, "Range, trc, vel, z, th, ph, sd, sa, ..." );

  pdef->addParam( "gaussian", "Width (in ms) at half maximum of gaussian to use instead of spike", NUM_VALUES_FIXED,
                  "0: spike on the nearest sample. A gaussian is evaluated at each sample" );
  pdef->addValue( "0", VALTYPE_NUMBER, "GAUSSIAN [ms]" );
  pdef->addParam( "tpower", "Hyperbola amplitudes scaled by (Time/t)**TPOWER", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "TPOWER" );
  pdef->addParam( "spj", "Shot point jump in units of GGD", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "SPJ" );
  pdef->addParam( "spi", "Shot point value of first record", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "SPI" );
  pdef->addParam( "zrngtrc", "Trace number of zero range trace", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "ZRNGTRC" );
  pdef->addParam( "vfile", "ASCII velocity file to convert to traces (one trace per record)", NUM_VALUES_FIXED,
                  "Lines: inline crossline X Y time velocity. A new trace starts when inline or crossline changes. HMR, HMT, GGD and the events are not used" );
  pdef->addValue( "NONE", VALTYPE_STRING, "File name" );
}

bool start_exec_mod_spk_synth_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}
void cleanup_mod_spk_synth_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  delete vars; vars = NULL;
}

extern "C" void _params_mod_spk_synth_( csParamDef* pdef ) {
  params_mod_spk_synth_( pdef );
}
extern "C" void _init_mod_spk_synth_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spk_synth_( param, env, writer );
}
extern "C" bool _start_exec_mod_spk_synth_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spk_synth_( env, writer );
}
extern "C" void _exec_mod_spk_synth_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_synth_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_synth_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spk_synth_( env, writer );
}
