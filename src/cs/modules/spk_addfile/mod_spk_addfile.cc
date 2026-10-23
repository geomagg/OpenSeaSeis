/* Copyright (c) 2026 M. A. Gallotti Guimaraes (COPPE-UFRJ).      */
/* Distributed under the OpenSeaSeis license (see LICENSE).       */
/* Module SPK_ADDFILE: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include "csSeismicReader.h"
#include <cmath>
#include <cstring>
#include <vector>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_ADDFILE
 *
 * SEISPAK process addfile (TYPE 2, existing OUTFILE = ONLINE1): adds the traces of
 * INFILE to the flow, sample by sample. Trace (record iy, trace ix) of the flow
 * receives trace (iy-IY0+1, ix-IX0+1) of INFILE; flow traces without a partner are
 * unchanged. With IX0 = IY0 = 1 (default) trace k of the flow gets trace k of INFILE.
 * For other IX0/IY0 the number of traces per record of the flow (nx_out, HMT) and of
 * INFILE (nx_in) must be given. INFILE samples beyond its length count as zero.
 * Not ported: creation of a new OUTFILE (NXOUT*NYOUT > 0) and TYPE 1/3 inputs
 * (transposed unformatted / dvsout files).
 */

namespace mod_spk_addfile {
  struct VariableStruct {
    csSeismicReader* reader; csTraceHeaderDef* hdef2; char* hdrBlock;
    int nt, nt1, ntrIn, ix0, iy0, nxOut, nxIn;
    long counter; int nadded;
    std::vector<float> buf;
    bool sequential;
  };
}
using namespace mod_spk_addfile;

void init_mod_spk_addfile_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* v = new VariableStruct();
  edef->setVariables( v );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  v->reader = NULL; v->hdef2 = NULL; v->hdrBlock = NULL;
  string fn;
  if( !param->exists("infile") ) writer->error("Parameter 'infile' is required");
  param->getString( "infile", &fn );
  v->ix0 = 1; v->iy0 = 1; v->nxOut = 0; v->nxIn = 0;
  if( param->exists("ix0") ) param->getInt( "ix0", &v->ix0 );
  if( param->exists("iy0") ) param->getInt( "iy0", &v->iy0 );
  if( param->exists("nx_out") ) param->getInt( "nx_out", &v->nxOut );
  if( param->exists("nx_in") )  param->getInt( "nx_in", &v->nxIn );
  v->sequential = ( v->ix0 == 1 && v->iy0 == 1 && ( v->nxOut <= 0 || v->nxOut == v->nxIn || v->nxIn <= 0 ) );
  if( !v->sequential && ( v->nxOut <= 0 || v->nxIn <= 0 ) )
    writer->error("IX0/IY0 other than 1 need nx_out (traces per record of the flow) and nx_in (of INFILE)");

  v->hdef2 = new csTraceHeaderDef( hdef );
  v->reader = new csSeismicReader( fn, !v->sequential, 1 );
  csSuperHeader shdr2; int blk = 0;
  if( !v->reader->readFileHeader( &shdr2, v->hdef2, &blk ) ) writer->error("Cannot read file header of '%s'", fn.c_str());
  v->hdef2->resetByteLocation();
  v->hdrBlock = new char[ blk + 8 ];
  v->nt = shdr->numSamples; v->nt1 = shdr2.numSamples; v->ntrIn = v->reader->numTraces();
  if( fabs( shdr2.sampleInt - shdr->sampleInt ) > 1e-4 )
    writer->warning("INFILE sample interval %g differs from the flow (%g): samples are added by index, as SEISPAK", shdr2.sampleInt, shdr->sampleInt);
  v->buf.assign( std::max( v->nt, v->nt1 ) + 8, 0.0f );
  v->counter = 0; v->nadded = 0;
  writer->line("  INFILE %s: %d traces, %d samples", fn.c_str(), v->ntrIn, v->nt1);
  if( v->sequential ) writer->line("  Trace k of the flow + trace k of INFILE");
  else writer->line("  IX0 %d IY0 %d, traces/record flow %d INFILE %d", v->ix0, v->iy0, v->nxOut, v->nxIn);
}

void exec_mod_spk_addfile_( csTraceGather* traceGather, int* port, int* numTrcToKeep,
                            csExecPhaseEnv* env, csLogWriter* writer )
{
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  float* s = traceGather->trace(0)->getTraceSamples();
  long k = v->counter++;
  long idx;
  if( v->sequential ) idx = k;
  else {
    int iy = (int)( k / v->nxOut ) + 1, ix = (int)( k % v->nxOut ) + 1;
    int irec1 = iy - v->iy0 + 1, itrc1 = ix - v->ix0 + 1;
    if( itrc1 < 1 || itrc1 > v->nxIn || irec1 < 1 ) return;
    idx = (long)(irec1-1)*v->nxIn + (itrc1-1);
  }
  if( idx >= v->ntrIn ) return;                                  // no partner: unchanged
  for( size_t j = 0; j < v->buf.size(); j++ ) v->buf[j] = 0.0f;
  if( !v->sequential ) v->reader->moveToTrace( (int)idx );
  if( !v->reader->readTrace( &v->buf[0], v->hdrBlock, v->nt1 ) ) return;
  for( int j = 0; j < v->nt; j++ ) s[j] += v->buf[j];
  v->nadded++;
}

void params_mod_spk_addfile_( csParamDef* pdef ) {
  pdef->setModule( "SPK_ADDFILE", "Add the traces of another file (SEISPAK addfile)",
                   "Adds INFILE trace by trace to the flow (record iy, trace ix gets INFILE trace iy-IY0+1, ix-IX0+1)." );
  pdef->addParam( "infile", "File to add (.cseis)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addParam( "ix0", "Start trace location in the flow", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "IX0" );
  pdef->addParam( "iy0", "Start record location in the flow", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "IY0" );
  pdef->addParam( "nx_out", "Traces per record of the flow (HMT), needed if IX0/IY0 != 1", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Traces" );
  pdef->addParam( "nx_in", "Traces per record of INFILE, needed if IX0/IY0 != 1", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "Traces" );
}
bool start_exec_mod_spk_addfile_( csExecPhaseEnv* env, csLogWriter* writer ) { return true; }
void cleanup_mod_spk_addfile_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( v == NULL ) return;
  writer->line("  SPK_ADDFILE: %d of %ld flow traces received an INFILE trace", v->nadded, v->counter);
  delete v->reader; delete v->hdef2; delete [] v->hdrBlock;
  delete v; v = NULL;
}
extern "C" void _params_mod_spk_addfile_( csParamDef* pdef ) { params_mod_spk_addfile_( pdef ); }
extern "C" void _init_mod_spk_addfile_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) { init_mod_spk_addfile_( param, env, writer ); }
extern "C" bool _start_exec_mod_spk_addfile_( csExecPhaseEnv* env, csLogWriter* writer ) { return start_exec_mod_spk_addfile_( env, writer ); }
extern "C" void _exec_mod_spk_addfile_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) { exec_mod_spk_addfile_( traceGather, port, numTrcToKeep, env, writer ); }
extern "C" void _cleanup_mod_spk_addfile_( csExecPhaseEnv* env, csLogWriter* writer ) { cleanup_mod_spk_addfile_( env, writer ); }
