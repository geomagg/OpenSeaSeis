/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module INPUT_SEISPAK: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <sys/types.h>

#include "csStandardHeaders.h"

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: INPUT_SEISPAK
 *
 * Reads seismic data in SEISPAK format (Texaco Seismic Process Development
 * Package, 'ice-pak' library, file signature 'ICEp').
 *
 * File layout (derived from ice-pak routines TGPSET, FX2FLT, GET1, GET2):
 *  - 12 byte preamble: 'ICEp', block size (words), number of blocks
 *  - Header directory: 4-byte words numbered from 1, word w at byte 12+4*(w-1)
 *     w1 IFW0 (first data word), w2 LHDR, w3 NT, w4 SMP [ms], w5 HMT (traces/record),
 *     w6 HMR (records), w9 NBYTE (0/4=float32, 2=scaled int16, 1=scaled byte),
 *     w19 ABIAS, w20 LNG2NRM (NBYTE=1 only), w21.. (name,value) pairs until 'END '
 *  - Traces, k = (ir-1)*HMT + (it-1):
 *     NBYTE 0/4: LHDR float headers + NT float32 samples
 *     NBYTE 2  : LHDR float headers + float FACTOR + [2 pad bytes if NT odd] + NT int16,
 *                sample = int16 / FACTOR
 *     NBYTE 1  : LHDR float headers + groups of LNG2NRM samples, each group:
 *                float FCTR + uint8 values, sample = (b-128)/FCTR + ABIAS
 */

namespace mod_input_seispak {

  static int const BASE_OFFSET = 12;

  struct SpkFile {
    std::string filename;
    FILE* fp;
    bool swap;
    int ifw0;
    int lhdr;
    int nt;
    float smp;
    int ihmt;
    int ihmr;
    int nbyte;
    float abias;
    int lng2nrm;
    long long recordBytes;
    long long firstTraceByte;
    int numTraces;
    std::vector<std::string> globalNames;
    std::vector<float> globalValues;
  };

  struct VariableStruct {
    int numFiles;
    SpkFile* files;
    int currentFile;
    int traceInFile;
    long traceCounter;
    long nTracesToRead;
    bool atEOF;
    char* buffer;
    int bufferSize;
    int numHdrOut;       // number of SeisPak trace headers exported
    int* hdrId_spk;      // SeaSeis header index for each SeisPak trace header
    int hdrId_trcno;
    int hdrId_fileno;
    int hdrId_rec;
    int hdrId_trc;
    int nSamplesOut;
  };

  static int const BYTEORDER_AUTO   = 0;
  static int const BYTEORDER_LITTLE = 1;
  static int const BYTEORDER_BIG    = 2;

  inline void swap4( char* p ) {
    char c;
    c = p[0]; p[0] = p[3]; p[3] = c;
    c = p[1]; p[1] = p[2]; p[2] = c;
  }
  inline void swap2( char* p ) {
    char c = p[0]; p[0] = p[1]; p[1] = c;
  }
  inline int getInt( char const* p, bool swap ) {
    char b[4]; memcpy( b, p, 4 ); if( swap ) swap4( b );
    int v; memcpy( &v, b, 4 ); return v;
  }
  inline float getFloat( char const* p, bool swap ) {
    char b[4]; memcpy( b, p, 4 ); if( swap ) swap4( b );
    float v; memcpy( &v, b, 4 ); return v;
  }
  inline short getShort( char const* p, bool swap ) {
    char b[2]; memcpy( b, p, 2 ); if( swap ) swap2( b );
    short v; memcpy( &v, b, 2 ); return v;
  }
  bool isLittleEndianMachine() {
    int one = 1;
    return( *(reinterpret_cast<char*>(&one)) == 1 );
  }

  // Open file, read header directory. Returns empty string on success, error message otherwise
  std::string openFile( SpkFile* f, int byteOrderOption ) {
    f->fp = fopen( f->filename.c_str(), "rb" );
    if( f->fp == NULL ) return "Cannot open file";
    char pre[BASE_OFFSET];
    if( fread( pre, 1, BASE_OFFSET, f->fp ) != (size_t)BASE_OFFSET ) return "File too short";
    if( strncmp( pre, "ICE", 3 ) != 0 ) return "SEISPAK signature 'ICE' not found at start of file";

    // First read words 1..20 to determine IFW0
    char dir[80];
    if( fread( dir, 1, 80, f->fp ) != 80 ) return "Header directory truncated";

    bool machineLittle = isLittleEndianMachine();
    if( byteOrderOption == BYTEORDER_AUTO ) {
      // LHDR (word 2) must be a small non-negative integer
      int lhdrNative = getInt( &dir[4], false );
      int lhdrSwap   = getInt( &dir[4], true );
      if( lhdrNative >= 0 && lhdrNative < 10000 ) f->swap = false;
      else if( lhdrSwap >= 0 && lhdrSwap < 10000 ) f->swap = true;
      else return "Cannot determine byte order (implausible LHDR value). Specify 'byte_order'";
    }
    else {
      bool fileLittle = ( byteOrderOption == BYTEORDER_LITTLE );
      f->swap = ( fileLittle != machineLittle );
    }
    bool sw = f->swap;
    f->ifw0  = getInt( &dir[0], sw );
    f->lhdr  = getInt( &dir[4], sw );
    f->nt    = getInt( &dir[8], sw );
    f->smp   = getFloat( &dir[12], sw );
    float hmt = getFloat( &dir[16], sw );
    float hmr = getFloat( &dir[20], sw );
    f->nbyte = getInt( &dir[32], sw );
    if( f->nbyte == 0 ) f->nbyte = 4;
    f->abias   = 0.0f;
    f->lng2nrm = 100;
    if( f->nbyte == 1 ) {
      f->abias = getFloat( &dir[72], sw );
      int l2n  = getInt( &dir[76], sw );
      if( l2n != 0 && l2n != 999 ) f->lng2nrm = l2n;
    }
    f->ihmt = (int)( hmt + 0.5f );
    f->ihmr = (int)( hmr + 0.5f );
    f->numTraces = f->ihmt * f->ihmr;

    if( f->ifw0 < 21 || f->lhdr < 0 || f->nt <= 0 || f->smp <= 0 || f->numTraces <= 0 ) {
      return "Implausible header directory values (IFW0/LHDR/NT/SMP/HMT/HMR)";
    }
    if( f->nbyte != 1 && f->nbyte != 2 && f->nbyte != 4 ) return "Unsupported NBYTE (word 9)";
    if( f->nbyte == 1 && f->lng2nrm <= 0 ) return "Invalid LNG2NRM (word 20)";

    // Global (name,value) pairs from word 21 up to IFW0-1
    int nwords = f->ifw0 - 1;
    std::vector<char> glob( 4*nwords );
    fseeko( f->fp, BASE_OFFSET, SEEK_SET );
    if( fread( &glob[0], 1, 4*nwords, f->fp ) != (size_t)(4*nwords) ) return "Header directory truncated";
    for( int w = 21; w+1 <= nwords; w += 2 ) {
      char const* p = &glob[4*(w-1)];
      if( strncmp( p, "END ", 4 ) == 0 ) break;
      std::string name;
      for( int i = 0; i < 4; i++ ) if( p[i] != ' ' && p[i] != 0 ) name += p[i];
      bool printable = !name.empty();
      for( size_t i = 0; i < name.size(); i++ ) if( name[i] < 33 || name[i] > 126 ) printable = false;
      if( !printable ) continue;
      f->globalNames.push_back( name );
      f->globalValues.push_back( getFloat( p+4, sw ) );
    }

    // Record size and position of first trace
    if( f->nbyte == 4 ) {
      f->recordBytes = 4LL * ( f->lhdr + f->nt );
    }
    else if( f->nbyte == 2 ) {
      int nnwd = f->lhdr + f->nt/2 + ( (f->nt % 2 == 0) ? 1 : 2 );
      f->recordBytes = 4LL * nnwd;
    }
    else {
      f->recordBytes = 4LL*f->lhdr + f->nt + 4LL*( (f->nt + f->lng2nrm - 1) / f->lng2nrm );
    }
    f->firstTraceByte = BASE_OFFSET + 4LL * ( f->ifw0 - 1 );
    return "";
  }

  // Read trace k (0-based) of file f into hdr[lhdr] and samples[nt]
  bool readTrace( SpkFile* f, int k, char* buffer, float* hdr, float* samples, int nSamplesOut ) {
    off_t pos = (off_t)( f->firstTraceByte + (long long)k * f->recordBytes );
    if( fseeko( f->fp, pos, SEEK_SET ) != 0 ) return false;
    if( fread( buffer, 1, (size_t)f->recordBytes, f->fp ) != (size_t)f->recordBytes ) return false;
    bool sw = f->swap;
    for( int i = 0; i < f->lhdr; i++ ) hdr[i] = getFloat( &buffer[4*i], sw );
    int nt = f->nt;
    int nCopy = ( nSamplesOut < nt ) ? nSamplesOut : nt;
    char const* p = &buffer[4*f->lhdr];

    if( f->nbyte == 4 ) {
      for( int i = 0; i < nCopy; i++ ) samples[i] = getFloat( p + 4*i, sw );
    }
    else if( f->nbyte == 2 ) {
      float factor = getFloat( p, sw );
      char const* s = p + 4 + ( (nt % 2) ? 2 : 0 );
      if( factor > 1.0e-30f ) {
        float inv = 1.0f / factor;
        for( int i = 0; i < nCopy; i++ ) samples[i] = inv * (float)getShort( s + 2*i, sw );
      }
      else {
        for( int i = 0; i < nCopy; i++ ) samples[i] = 0.0f;  // same behaviour as GET2
      }
    }
    else { // nbyte == 1
      int l2n = f->lng2nrm;
      int isamp = 0;
      char const* g = p;
      while( isamp < nt ) {
        int n = ( nt - isamp < l2n ) ? ( nt - isamp ) : l2n;
        float fctr = getFloat( g, sw );
        float inv = ( fctr != 0.0f ) ? 1.0f/fctr : 0.0f;
        unsigned char const* b = reinterpret_cast<unsigned char const*>( g + 4 );
        for( int j = 0; j < n; j++ ) {
          if( isamp + j < nCopy ) samples[isamp+j] = ( (float)b[j] - 128.0f ) * inv + f->abias;
        }
        isamp += n;
        g += 4 + n;
      }
    }
    for( int i = nCopy; i < nSamplesOut; i++ ) samples[i] = 0.0f;
    return true;
  }

  std::string toLower( std::string const& s ) {
    std::string r = s;
    for( size_t i = 0; i < r.size(); i++ ) r[i] = (char)tolower( r[i] );
    return r;
  }
}

using namespace mod_input_seispak;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_input_seispak_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_INPUT );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->numFiles      = 0;
  vars->files         = NULL;
  vars->currentFile   = 0;
  vars->traceInFile   = 0;
  vars->traceCounter  = 0;
  vars->nTracesToRead = -1;
  vars->atEOF         = false;
  vars->buffer        = NULL;
  vars->bufferSize    = 0;
  vars->numHdrOut     = 0;
  vars->hdrId_spk     = NULL;
  vars->hdrId_trcno   = -1;
  vars->hdrId_fileno  = -1;
  vars->hdrId_rec     = -1;
  vars->hdrId_trc     = -1;
  vars->nSamplesOut   = 0;

  std::string text;

  //---------------------------------------------
  int byteOrder = BYTEORDER_AUTO;
  if( param->exists("byte_order") ) {
    param->getString( "byte_order", &text );
    if( !text.compare("auto") )        byteOrder = BYTEORDER_AUTO;
    else if( !text.compare("little") ) byteOrder = BYTEORDER_LITTLE;
    else if( !text.compare("big") )    byteOrder = BYTEORDER_BIG;
    else writer->error("Unknown option for 'byte_order': %s", text.c_str());
  }

  if( param->exists("ntraces") ) {
    int n = 0;
    param->getInt( "ntraces", &n );
    vars->nTracesToRead = ( n > 0 ) ? n : -1;
  }

  bool useGlobalNames = true;
  if( param->exists("hdr_names") ) {
    param->getString( "hdr_names", &text );
    if( !text.compare("yes") ) useGlobalNames = true;
    else if( !text.compare("no") ) useGlobalNames = false;
    else writer->error("Unknown option for 'hdr_names': %s", text.c_str());
  }

  //---------------------------------------------
  vars->numFiles = param->getNumLines("filename");
  if( vars->numFiles <= 0 ) writer->error("Required input parameter missing: 'filename'");
  vars->files = new SpkFile[vars->numFiles];
  for( int ifile = 0; ifile < vars->numFiles; ifile++ ) {
    SpkFile* f = &vars->files[ifile];
    f->fp = NULL;
    param->getStringAtLine( "filename", &f->filename, ifile );
    std::string err = openFile( f, byteOrder );
    if( !err.empty() ) writer->error("Error reading SEISPAK file '%s': %s", f->filename.c_str(), err.c_str());
    if( ifile > 0 ) {
      SpkFile* f0 = &vars->files[0];
      if( f->nt != f0->nt ) writer->error("File '%s': number of samples (%d) differs from first file (%d)", f->filename.c_str(), f->nt, f0->nt);
      if( f->smp != f0->smp ) writer->error("File '%s': sample interval (%f) differs from first file (%f)", f->filename.c_str(), f->smp, f0->smp);
      if( f->lhdr != f0->lhdr ) writer->error("File '%s': number of trace headers (%d) differs from first file (%d)", f->filename.c_str(), f->lhdr, f0->lhdr);
    }
    if( (int)f->recordBytes > vars->bufferSize ) vars->bufferSize = (int)f->recordBytes;
  }
  vars->buffer = new char[vars->bufferSize];

  SpkFile* f0 = &vars->files[0];
  int nsamples = 0;
  if( param->exists("nsamples") ) param->getInt( "nsamples", &nsamples );
  vars->nSamplesOut = ( nsamples > 0 ) ? nsamples : f0->nt;
  shdr->numSamples = vars->nSamplesOut;
  shdr->sampleInt  = f0->smp;

  //---------------------------------------------
  // Log
  char const* fmtName[5] = { "", "scaled byte", "scaled int16", "", "float32" };
  writer->line("");
  for( int ifile = 0; ifile < vars->numFiles; ifile++ ) {
    SpkFile* f = &vars->files[ifile];
    writer->line("  File %2d: %s", ifile+1, f->filename.c_str());
    writer->line("     Byte order swapped: %s,  NBYTE: %d (%s)", f->swap ? "yes" : "no", f->nbyte, fmtName[f->nbyte]);
    writer->line("     IFW0: %d  LHDR: %d  NT: %d  SMP: %f ms", f->ifw0, f->lhdr, f->nt, f->smp);
    writer->line("     Traces/record (HMT): %d  Records (HMR): %d  Total traces: %d", f->ihmt, f->ihmr, f->numTraces);
    for( size_t i = 0; i < f->globalNames.size(); i++ ) {
      writer->line("     Global: %-4s = %g", f->globalNames[i].c_str(), f->globalValues[i]);
    }
  }
  writer->line("");

  //---------------------------------------------
  // Trace headers
  vars->hdrId_trcno  = hdef->addStandardHeader( HDR_TRCNO.name );
  vars->hdrId_fileno = hdef->addStandardHeader( HDR_FILENO.name );
  vars->hdrId_rec    = hdef->addHeader( TYPE_INT, "spk_rec", "SEISPAK record number (1..HMR)" );
  vars->hdrId_trc    = hdef->addHeader( TYPE_INT, "spk_trc", "SEISPAK trace number within record (1..HMT)" );

  vars->numHdrOut = f0->lhdr;
  vars->hdrId_spk = new int[vars->numHdrOut > 0 ? vars->numHdrOut : 1];
  for( int ih = 0; ih < vars->numHdrOut; ih++ ) {
    char defName[32];
    sprintf( defName, "spk_h%d", ih+1 );
    std::string name = defName;
    std::string desc = "SEISPAK trace header word " + std::string(defName+5);
    if( useGlobalNames ) {
      // Global pair (NAME, position) assigns a name to trace header word 'position'
      for( size_t ig = 0; ig < f0->globalNames.size(); ig++ ) {
        float v = f0->globalValues[ig];
        if( fabs( v - (float)(ih+1) ) < 1.0e-4 ) {
          name = "spk_" + toLower( f0->globalNames[ig] );
          desc = "SEISPAK trace header " + f0->globalNames[ig];
          break;
        }
      }
    }
    if( hdef->headerExists( name ) ) name = std::string(defName);
    vars->hdrId_spk[ih] = hdef->addHeader( TYPE_FLOAT, name, desc );
    writer->line("  Trace header word %2d  ->  SeaSeis header '%s'", ih+1, name.c_str());
  }
  writer->line("");
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_input_seispak_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );

  if( vars->atEOF || ( vars->nTracesToRead > 0 && vars->traceCounter >= vars->nTracesToRead ) ) {
    vars->atEOF = true;
    traceGather->freeAllTraces();
    return;
  }

  // Advance to next file if current one is exhausted
  while( vars->currentFile < vars->numFiles && vars->traceInFile >= vars->files[vars->currentFile].numTraces ) {
    SpkFile* f = &vars->files[vars->currentFile];
    if( f->fp != NULL ) { fclose( f->fp ); f->fp = NULL; }
    vars->currentFile += 1;
    vars->traceInFile = 0;
  }
  if( vars->currentFile >= vars->numFiles ) {
    vars->atEOF = true;
    traceGather->freeAllTraces();
    return;
  }

  SpkFile* f = &vars->files[vars->currentFile];
  csTrace* trace = traceGather->trace(0);
  csTraceHeader* trcHdr = trace->getTraceHeader();
  float* samples = trace->getTraceSamples();

  std::vector<float> hdr( f->lhdr > 0 ? f->lhdr : 1 );
  int k = vars->traceInFile;
  if( !readTrace( f, k, vars->buffer, &hdr[0], samples, vars->nSamplesOut ) ) {
    writer->error("Error reading trace %d from SEISPAK file '%s' (file truncated?)", k+1, f->filename.c_str());
  }

  trcHdr->setIntValue( vars->hdrId_trcno, (int)vars->traceCounter + 1 );
  trcHdr->setIntValue( vars->hdrId_fileno, vars->currentFile + 1 );
  trcHdr->setIntValue( vars->hdrId_rec, k / f->ihmt + 1 );
  trcHdr->setIntValue( vars->hdrId_trc, k % f->ihmt + 1 );
  for( int ih = 0; ih < vars->numHdrOut; ih++ ) {
    trcHdr->setFloatValue( vars->hdrId_spk[ih], hdr[ih] );
  }

  vars->traceInFile  += 1;
  vars->traceCounter += 1;
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_input_seispak_( csParamDef* pdef ) {
  pdef->setModule( "INPUT_SEISPAK", "Input data in SEISPAK format (Texaco ice-pak, 'ICEp' files)",
                   "Reads SEISPAK disk files written by the ice-pak library. Supported sample formats (header word 9, NBYTE): "
                   "0/4 = 32bit float, 2 = 16bit integer scaled per trace, 1 = 8bit integer scaled per sample group. "
                   "Traces are output in file order (record by record). "
                   "SEISPAK trace header words are stored in float headers 'spk_h1', 'spk_h2'... "
                   "or 'spk_<name>' if the file's global header list assigns a name to that word position (e.g. RNUM, TNUM). "
                   "Record and trace numbers are stored in 'spk_rec' and 'spk_trc'." );

  pdef->addParam( "filename", "Input file name", NUM_VALUES_FIXED, "Specify one line per input file. Files are read in sequence." );
  pdef->addValue( "", VALTYPE_STRING, "Input file name" );

  pdef->addParam( "nsamples", "Number of samples to read in", NUM_VALUES_FIXED,
                  "Traces are truncated or padded with zeros. Set 0 to use number of samples of input file.");
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of samples" );

  pdef->addParam( "ntraces", "Number of traces to read in", NUM_VALUES_FIXED, "Set 0 to read all traces." );
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of traces" );

  pdef->addParam( "byte_order", "Byte order of input file", NUM_VALUES_FIXED );
  pdef->addValue( "auto", VALTYPE_OPTION );
  pdef->addOption( "auto", "Detect byte order from header directory" );
  pdef->addOption( "little", "Little endian (e.g. DEC, Linux PC)" );
  pdef->addOption( "big", "Big endian (e.g. SGI, SUN, HP, IBM)" );

  pdef->addParam( "hdr_names", "Use global header names for trace header words", NUM_VALUES_FIXED );
  pdef->addValue( "yes", VALTYPE_OPTION );
  pdef->addOption( "yes", "Name trace headers 'spk_<name>' from global (name,position) list" );
  pdef->addOption( "no", "Name trace headers 'spk_h1', 'spk_h2', ..." );
}

//************************************************************************************************
// Start exec phase
//*************************************************************************************************
bool start_exec_mod_input_seispak_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

//************************************************************************************************
// Cleanup phase
//*************************************************************************************************
void cleanup_mod_input_seispak_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( vars->files != NULL ) {
    for( int ifile = 0; ifile < vars->numFiles; ifile++ ) {
      if( vars->files[ifile].fp != NULL ) {
        fclose( vars->files[ifile].fp );
        vars->files[ifile].fp = NULL;
      }
    }
    delete [] vars->files;
    vars->files = NULL;
  }
  if( vars->buffer != NULL ) { delete [] vars->buffer; vars->buffer = NULL; }
  if( vars->hdrId_spk != NULL ) { delete [] vars->hdrId_spk; vars->hdrId_spk = NULL; }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_input_seispak_( csParamDef* pdef ) {
  params_mod_input_seispak_( pdef );
}
extern "C" void _init_mod_input_seispak_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_input_seispak_( param, env, writer );
}
extern "C" bool _start_exec_mod_input_seispak_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_input_seispak_( env, writer );
}
extern "C" void _exec_mod_input_seispak_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_input_seispak_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_input_seispak_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_input_seispak_( env, writer );
}
