/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module SPK_VXT2ONL: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

#include "csStandardHeaders.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_VXT2ONL
 *
 * Input module: converts a SEISPAK VXT velocity list (horizons, as in the DVSMIG
 * jobdeck) into interval-velocity traces v(z), one trace per DPN, sampled in depth.
 * Write them with OUTPUT_SEISPAK (format float) to obtain an ONL velocity file
 * (as used by DVSMIG VEL -3 / ONLVELS and by SPK_DVSMIG2D 'vel_file').
 * The VXT list is converted by the original VEVENTS routine (shared with SPK_DVSMIG2D).
 */

extern "C" {
  void spkvcol_( float* xx, int* np, int* nz, int* itzr, float* ddz, float* dt, int* irfc, float* e );
  void spkvxt_( float* vel, float* dx, int* nxv, int* nyv, int* np, float* afit, float* ee, float* x, float* bfit, float* xx );
}

namespace mod_spk_vxt2onl {
  struct VariableStruct {
    std::vector<float> afit;
    int nxv, np;
    int nx;            // number of traces to output
    int nnz;
    float ddz;
    int traceCounter;
    bool atEOF;
    int hdrId_cdp;
    int hdrId_trcno;
  };
}
using mod_spk_vxt2onl::VariableStruct;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_spk_vxt2onl_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  edef->setExecType( EXEC_TYPE_INPUT );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->traceCounter = 0;
  vars->atEOF = false;
  vars->nx = 0;
  vars->nnz = 0;
  vars->ddz = 0.0f;

  std::string fname;
  param->getString( "vxt_file", &fname );
  param->getFloat( "ddz", &vars->ddz );
  param->getInt( "nnz", &vars->nnz );
  if( vars->ddz <= 0.0f || vars->nnz <= 0 ) writer->error("DDZ and NNZ must be > 0");
  if( param->exists("nx") ) param->getInt( "nx", &vars->nx );
  float dx = 1.0f;
  if( param->exists("dx") ) param->getFloat( "dx", &dx );
  float yLine = -1.0f;
  if( param->exists("vxt_line") ) param->getFloat( "vxt_line", &yLine );

  //--- read numbers after the keyword VXT (or VEL); if no keyword, the whole file
  FILE* fp = fopen( fname.c_str(), "r" );
  if( fp == NULL ) writer->error("Cannot open VXT file '%s'", fname.c_str());
  std::vector<std::string> tok;
  char word[256];
  while( fscanf( fp, "%255s", word ) == 1 ) tok.push_back( word );
  fclose( fp );
  size_t start = 0;
  for( size_t i = 0; i < tok.size(); i++ ) {
    std::string t = tok[i];
    for( size_t c = 0; c < t.size(); c++ ) t[c] = (char)tolower(t[c]);
    if( t == "vxt" || t == "vel" ) { start = i+1; break; }
  }
  std::vector<float> v;
  for( size_t i = start; i < tok.size(); i++ ) {
    char* end = NULL;
    float f = strtof( tok[i].c_str(), &end );
    if( end == tok[i].c_str() || *end != 0 ) break;
    v.push_back( f );
  }
  if( v.size() < 5 ) writer->error("VXT file '%s': no velocity list found", fname.c_str());

  //--- split into y-line groups: [code] [negative limits] y1 <horizons...> 0 y2 <horizons...> 0 ...
  size_t pos = 1;
  std::vector<float> lead;
  while( pos < v.size() && v[pos] < 0.0f ) lead.push_back( v[pos++] );
  std::vector<float> ys;
  std::vector< std::vector<float> > bodies;
  while( pos < v.size() ) {
    float y = v[pos++];
    size_t s0 = pos;
    while( pos+1 < v.size() && !( v[pos] <= 0.5f && v[pos+1] <= 0.5f ) ) pos++;
    bodies.push_back( std::vector<float>( v.begin()+s0, v.begin()+std::min(pos+1, v.size()) ) );
    ys.push_back( y );
    pos += 2;
  }
  int ig = 0;
  if( yLine >= 0.0f ) {
    ig = -1;
    for( size_t i = 0; i < ys.size(); i++ ) if( fabs( ys[i] - yLine ) < 0.01f ) ig = (int)i;
    if( ig < 0 ) writer->error("VXT file: line y=%g not found (%d lines, y = %g ... %g)", yLine, (int)ys.size(), ys.front(), ys.back());
  }

  //--- 2D list for VEVENTS: code, limits, 1, horizons of the selected line, 0
  std::vector<float> list;
  list.push_back( v[0] );
  list.insert( list.end(), lead.begin(), lead.end() );
  list.push_back( 1.0f );
  list.insert( list.end(), bodies[ig].begin(), bodies[ig].end() );
  list.push_back( 0.0f );
  size_t m = list.size() + 60000;
  list.resize( m, 0.0f );
  std::vector<float> ee( m ), x( 20*m + 100000 ), xx( 10001 );
  size_t nafit = 3 + (size_t)10001*203;   // up to 10000 DPN and 100 horizons
  vars->afit.assign( nafit, 0.0f );
  std::vector<float> bfit( nafit );
  int nxv = 0, nyv = 0, np = 0;
  spkvxt_( &list[0], &dx, &nxv, &nyv, &np, &vars->afit[0], &ee[0], &x[0], &bfit[0], &xx[0] );
  if( nxv <= 0 || np <= 0 || nxv > 10000 ) writer->error("VXT file '%s': could not build velocity functions (nx=%d, events=%d)", fname.c_str(), nxv, np);
  vars->nxv = nxv;
  vars->np  = np;
  if( vars->nx <= 0 ) vars->nx = nxv;

  shdr->numSamples = vars->nnz;
  shdr->sampleInt  = vars->ddz;
  shdr->domain     = DOMAIN_XD;

  vars->hdrId_trcno = hdef->addStandardHeader( HDR_TRCNO.name );
  vars->hdrId_cdp   = hdef->addStandardHeader( HDR_CMP.name );

  writer->line("  VXT file: %s  (%d y-lines in file, using y = %g)", fname.c_str(), (int)ys.size(), ys[ig]);
  writer->line("  VXT code %g, %d horizons, %d locations (DPN)", v[0], np, nxv);
  writer->line("  Output: %d traces (DPN 1..%d), %d depth samples, DDZ %f (max. depth %f)",
               vars->nx, vars->nx, vars->nnz, vars->ddz, (vars->nnz-1)*vars->ddz);
  writer->line("  Sample j = interval velocity between depths (j-1)*DDZ and j*DDZ (as onl2dvsvels)");
  if( vars->nx > nxv ) writer->warning("NX (%d) larger than the VXT extent (%d DPN): the last velocity function is repeated", vars->nx, nxv);
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_spk_vxt2onl_(
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
  int k = vars->traceCounter;
  if( k == vars->nx-1 ) vars->atEOF = true;

  csTrace* trace = traceGather->trace(0);
  float* samples = trace->getTraceSamples();
  int ix = k + 1;
  if( ix > vars->nxv ) ix = vars->nxv;
  int np = vars->np, nw = 2*np+3;
  int nnz = vars->nnz, itzr = 1, irfc = 1;
  float dt = 0.0f;
  float* xx = &vars->afit[3 + (size_t)(ix-1)*nw];
  spkvcol_( xx, &np, &nnz, &itzr, &vars->ddz, &dt, &irfc, samples );

  csTraceHeader* trcHdr = trace->getTraceHeader();
  trcHdr->setIntValue( vars->hdrId_trcno, k+1 );
  trcHdr->setIntValue( vars->hdrId_cdp, k+1 );
  vars->traceCounter += 1;
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_spk_vxt2onl_( csParamDef* pdef ) {
  pdef->setModule( "SPK_VXT2ONL", "Convert SEISPAK VXT velocity list to v(z) traces (ONL velocity file)",
                   "Input module. Generates one interval-velocity trace per DPN, sampled in depth at DDZ, from a SEISPAK VXT "
                   "horizon list (as in the DVSMIG jobdeck), using the original VEVENTS routine. "
                   "Use OUTPUT_SEISPAK with 'format float' to write the ONL file. Header 'cmp' = DPN." );

  pdef->addParam( "vxt_file", "File with SEISPAK VXT velocity list", NUM_VALUES_FIXED,
                  "Numbers after the keyword VXT (or the whole file): code, then for each y-line: y, horizons as triplets "
                  "(velocity, DPN, depth or time) each ended by a value <= 0.5, and 0 at the end of the line." );
  pdef->addValue( "", VALTYPE_STRING, "File name" );

  pdef->addParam( "vxt_line", "y-line of the VXT file to use (3D models)", NUM_VALUES_FIXED, "Default: first line in the file" );
  pdef->addValue( "", VALTYPE_NUMBER, "y value as written in the file" );

  pdef->addParam( "ddz", "Depth sample rate of the output", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "DDZ" );

  pdef->addParam( "nnz", "Number of depth samples", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NNZ" );

  pdef->addParam( "nx", "Number of output traces (DPN 1..NX)", NUM_VALUES_FIXED, "Default: DPN extent of the VXT list" );
  pdef->addValue( "0", VALTYPE_NUMBER, "NX" );

  pdef->addParam( "dx", "Distance between DPNs (only for VXT codes >= 5, positions given in distance)", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "DX" );
}

//************************************************************************************************
bool start_exec_mod_spk_vxt2onl_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

void cleanup_mod_spk_vxt2onl_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  delete vars; vars = NULL;
}

extern "C" void _params_mod_spk_vxt2onl_( csParamDef* pdef ) {
  params_mod_spk_vxt2onl_( pdef );
}
extern "C" void _init_mod_spk_vxt2onl_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spk_vxt2onl_( param, env, writer );
}
extern "C" bool _start_exec_mod_spk_vxt2onl_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spk_vxt2onl_( env, writer );
}
extern "C" void _exec_mod_spk_vxt2onl_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_vxt2onl_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_vxt2onl_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spk_vxt2onl_( env, writer );
}
