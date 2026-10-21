/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module SPK_MFILT2: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstring>
#include <cstdio>
#include <vector>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_MFILT2
 *
 * SEISPAK program mfilt2: reads the filter in FILTER (Fortran unformatted:
 * record 1 = NF NTRF NRF SRF, record 2 = NF coefficients, as written by lsqmatch /
 * SPK_LSQMATCH) and convolves every trace with it, centred on sample NF/2+1
 * (CONVG, output samples NF/2+1 .. NF/2+NT).
 */

namespace mod_spk_mfilt2 {
  struct VariableStruct { std::vector<float> f; int nf; };
}
using namespace mod_spk_mfilt2;

static bool read_rec( FILE* fi, std::vector<char>& buf ) {
  int n1 = 0, n2 = 0;
  if( fread( &n1, 4, 1, fi ) != 1 || n1 < 0 || n1 > 100000000 ) return false;
  buf.resize( n1 );
  if( n1 > 0 && fread( &buf[0], 1, n1, fi ) != (size_t)n1 ) return false;
  if( fread( &n2, 4, 1, fi ) != 1 || n2 != n1 ) return false;
  return true;
}

void init_mod_spk_mfilt2_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef* edef = env->execPhaseDef;
  VariableStruct* v = new VariableStruct();
  edef->setVariables( v );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  string fn = "FILTER";
  if( param->exists("filter_file") ) param->getString( "filter_file", &fn );
  FILE* fi = fopen( fn.c_str(), "rb" );
  if( !fi ) writer->error("Cannot open filter file '%s'", fn.c_str());
  std::vector<char> b;
  if( !read_rec( fi, b ) || b.size() < 16 ) writer->error("'%s': first record (NF NTRF NRF SRF) not readable", fn.c_str());
  int nf, ntrf, nrf; float srf;
  memcpy( &nf, &b[0], 4 ); memcpy( &ntrf, &b[4], 4 ); memcpy( &nrf, &b[8], 4 ); memcpy( &srf, &b[12], 4 );
  if( nf < 1 || nf > 100000 ) writer->error("'%s': invalid filter length %d", fn.c_str(), nf);
  if( !read_rec( fi, b ) || (int)b.size() < 4*nf ) writer->error("'%s': filter record not readable", fn.c_str());
  fclose( fi );
  v->nf = nf; v->f.resize( nf );
  memcpy( &v->f[0], &b[0], 4*nf );
  writer->line("  filter %s: nf %d ntrf %d nrf %d srf %g", fn.c_str(), nf, ntrf, nrf, srf);
  writer->line("  centred on coefficient %d", nf/2+1);
}

void exec_mod_spk_mfilt2_( csTraceGather* traceGather, int* port, int* numTrcToKeep,
                           csExecPhaseEnv* env, csLogWriter* writer )
{
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  int nt = env->superHeader->numSamples;
  float* s = traceGather->trace(0)->getTraceSamples();
  int nf = v->nf, nf2 = nf/2;
  std::vector<float> out( nt, 0.0f );
  // CONVG(f, nf, tr, nt, c, LA = nf2+1, LB = nf2+nt+1): out(k) = c(nf2+k)
  for( int k = 1; k <= nt; k++ ) {
    int jj = nf2 + k;                                // conv index (1-based): sum f(K) tr(I), K+I-1 = jj
    double acc = 0.0;
    for( int i = std::max( 1, jj-nf+1 ); i <= std::min( nt, jj ); i++ ) acc += v->f[jj-i]*s[i-1];
    out[k-1] = (float)acc;
  }
  memcpy( s, &out[0], nt*sizeof(float) );
}

void params_mod_spk_mfilt2_( csParamDef* pdef ) {
  pdef->setModule( "SPK_MFILT2", "Apply the filter in FILTER (SEISPAK mfilt2)",
                   "Convolves every trace with the filter written by SPK_LSQMATCH (or SEISPAK lsqmatch), centred on coefficient NF/2+1." );
  pdef->addParam( "filter_file", "Filter file (Fortran unformatted, SEISPAK FILTER format)", NUM_VALUES_FIXED );
  pdef->addValue( "FILTER", VALTYPE_STRING, "File name" );
}
bool start_exec_mod_spk_mfilt2_( csExecPhaseEnv* env, csLogWriter* writer ) { return true; }
void cleanup_mod_spk_mfilt2_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  delete v; v = NULL;
}
extern "C" void _params_mod_spk_mfilt2_( csParamDef* pdef ) { params_mod_spk_mfilt2_( pdef ); }
extern "C" void _init_mod_spk_mfilt2_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) { init_mod_spk_mfilt2_( param, env, writer ); }
extern "C" bool _start_exec_mod_spk_mfilt2_( csExecPhaseEnv* env, csLogWriter* writer ) { return start_exec_mod_spk_mfilt2_( env, writer ); }
extern "C" void _exec_mod_spk_mfilt2_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) { exec_mod_spk_mfilt2_( traceGather, port, numTrcToKeep, env, writer ); }
extern "C" void _cleanup_mod_spk_mfilt2_( csExecPhaseEnv* env, csLogWriter* writer ) { cleanup_mod_spk_mfilt2_( env, writer ); }
