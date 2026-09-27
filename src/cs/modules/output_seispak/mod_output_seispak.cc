/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module OUTPUT_SEISPAK: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <sys/types.h>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: OUTPUT_SEISPAK
 *
 * Writes seismic data in SEISPAK format (Texaco Seismic Process Development
 * Package, 'ice-pak' library, file signature 'ICEp').
 * Encoding follows ice-pak routines TGPSET, FLT2FX, PUT1 and PUT2.
 * See module INPUT_SEISPAK for a description of the file layout.
 *
 * The header directory is written at the end (cleanup phase), once the total
 * number of traces is known.
 */

namespace mod_output_seispak {

  static int const BASE_OFFSET = 12;
  static int const BLOCK_WORDS = 1024;        // ice-pak block size in 4-byte words
  static int const NULL_VALUE  = 999;         // ice-pak 'undefined' value
  static float const EPSLON    = 1.0e-30f;

  static int const SRC_HEADER = 0;   // value taken from SeaSeis trace header
  static int const SRC_RECORD = 1;   // SEISPAK record number
  static int const SRC_TRACE  = 2;   // SEISPAK trace number within record

  struct HdrOut {
    int source;
    int hdrId;
    std::string spkName;   // up to 4 characters; empty = unnamed word
  };

  struct VariableStruct {
    std::string filename;
    FILE* fp;
    bool swap;
    int nbyte;
    int lng2nrm;
    int ntracesPerRecord;   // 0 = all traces in one record
    int nt;
    float smp;
    std::vector<HdrOut> hdrs;
    int ifw0;
    long long recordBytes;
    long long firstTraceByte;
    long traceCounter;
    char* buffer;
    std::vector<float> work;
    bool isOK;
  };

  bool isLittleEndianMachine() {
    int one = 1;
    return( *(reinterpret_cast<char*>(&one)) == 1 );
  }
  inline void putInt( char* p, int v, bool swap ) {
    memcpy( p, &v, 4 );
    if( swap ) { char c; c=p[0]; p[0]=p[3]; p[3]=c; c=p[1]; p[1]=p[2]; p[2]=c; }
  }
  inline void putFloat( char* p, float v, bool swap ) {
    memcpy( p, &v, 4 );
    if( swap ) { char c; c=p[0]; p[0]=p[3]; p[3]=c; c=p[1]; p[1]=p[2]; p[2]=c; }
  }
  inline void putShort( char* p, short v, bool swap ) {
    memcpy( p, &v, 2 );
    if( swap ) { char c=p[0]; p[0]=p[1]; p[1]=c; }
  }
  std::string toUpper4( std::string const& s ) {
    std::string r;
    for( size_t i = 0; i < s.size() && r.size() < 4; i++ ) r += (char)toupper( s[i] );
    return r;
  }

  // Encode one trace (headers + samples) into buffer, following PUT2 / PUT1 / FLT2FX
  void encodeTrace( VariableStruct* vars, float const* hdr, float const* samples, char* buf ) {
    int lhdr = (int)vars->hdrs.size();
    int nt   = vars->nt;
    bool sw  = vars->swap;
    for( int i = 0; i < lhdr; i++ ) putFloat( &buf[4*i], hdr[i], sw );
    char* p = &buf[4*lhdr];

    if( vars->nbyte == 4 ) {
      for( int i = 0; i < nt; i++ ) putFloat( p + 4*i, samples[i], sw );
      return;
    }
    if( vars->nbyte == 2 ) {
      float amax = 0.0f;
      for( int i = 0; i < nt; i++ ) { float a = fabsf( samples[i] ); if( a > amax ) amax = a; }
      if( amax < EPSLON || amax != amax ) amax = 1.0f;
      float factor = 32000.0f / amax;
      putFloat( p, factor, sw );
      char* s = p + 4;
      if( nt % 2 ) { putShort( s, 0, sw ); s += 2; }   // pad value (undefined in ice-pak)
      // Rounded (ice-pak truncates): avoids 1-LSB loss when re-writing int16 data read from SEISPAK
      for( int i = 0; i < nt; i++ ) putShort( s + 2*i, (short)lrintf( factor * samples[i] ), sw );
      return;
    }
    // nbyte == 1
    int l2n = vars->lng2nrm;
    char* g = p;
    for( int i0 = 0; i0 < nt; i0 += l2n ) {
      int n = ( nt - i0 < l2n ) ? ( nt - i0 ) : l2n;
      float amx = 0.0f;
      for( int j = 0; j < n; j++ ) { float a = fabsf( samples[i0+j] ); if( a > amx ) amx = a; }
      float fctr = ( amx > 0.0f ) ? 127.49f / amx : 1.0f;
      putFloat( g, fctr, sw );
      unsigned char* b = reinterpret_cast<unsigned char*>( g + 4 );
      for( int j = 0; j < n; j++ ) {
        int iv = (int)( 128.5f + samples[i0+j] * fctr );
        if( iv < 0 ) iv = 0;
        if( iv > 255 ) iv = 255;
        b[j] = (unsigned char)iv;
      }
      g += 4 + n;
    }
  }

  // Write header directory + preamble and pad file to a full block. Called at end.
  bool finalizeFile( VariableStruct* vars, csLogWriter* writer ) {
    if( vars->fp == NULL ) return false;
    long ntr = vars->traceCounter;
    int hmt = ( vars->ntracesPerRecord > 0 ) ? vars->ntracesPerRecord : (int)ntr;
    if( hmt <= 0 ) hmt = 1;
    int hmr = (int)( ( ntr + hmt - 1 ) / hmt );
    if( hmr <= 0 ) hmr = 1;
    int lhdr = (int)vars->hdrs.size();

    // Fill last record with zero traces
    long ntrTotal = (long)hmt * hmr;
    if( ntrTotal > ntr ) {
      std::vector<float> hdr( lhdr > 0 ? lhdr : 1, 0.0f );
      std::vector<float> zero( vars->nt, 0.0f );
      for( long k = ntr; k < ntrTotal; k++ ) {
        for( int ih = 0; ih < lhdr; ih++ ) {
          if( vars->hdrs[ih].source == SRC_RECORD ) hdr[ih] = (float)( k / hmt + 1 );
          else if( vars->hdrs[ih].source == SRC_TRACE ) hdr[ih] = (float)( k % hmt + 1 );
          else hdr[ih] = 0.0f;
        }
        encodeTrace( vars, &hdr[0], &zero[0], vars->buffer );
        off_t pos = (off_t)( vars->firstTraceByte + (long long)k * vars->recordBytes );
        fseeko( vars->fp, pos, SEEK_SET );
        fwrite( vars->buffer, 1, (size_t)vars->recordBytes, vars->fp );
      }
      writer->line("  Last record padded with %ld zero trace(s)", ntrTotal - ntr);
    }

    // Pad to full block with NULL_VALUE words
    long long endByte  = vars->firstTraceByte + (long long)ntrTotal * vars->recordBytes;
    long long blkBytes = 4LL * BLOCK_WORDS;
    long long nblk     = ( endByte + blkBytes - 1 ) / blkBytes;
    long long padBytes = nblk * blkBytes - endByte;
    fseeko( vars->fp, (off_t)endByte, SEEK_SET );
    char w999[4];
    putInt( w999, NULL_VALUE, vars->swap );
    for( long long i = 0; i < padBytes / 4; i++ ) fwrite( w999, 1, 4, vars->fp );
    for( long long i = 0; i < padBytes % 4; i++ ) fputc( 0, vars->fp );

    // Preamble + header directory
    int nwords = vars->ifw0 - 1;
    std::vector<char> dir( BASE_OFFSET + 4*nwords );
    memcpy( &dir[0], "ICEp", 4 );
    putInt( &dir[4], BLOCK_WORDS, vars->swap );
    putInt( &dir[8], (int)nblk, vars->swap );
    char* d = &dir[BASE_OFFSET];
    for( int w = 1; w <= nwords; w++ ) putInt( &d[4*(w-1)], NULL_VALUE, vars->swap );
    putInt(   &d[0],  vars->ifw0, vars->swap );     // w1 IFW0
    putInt(   &d[4],  lhdr,       vars->swap );     // w2 LHDR
    putInt(   &d[8],  vars->nt,   vars->swap );     // w3 NT
    putFloat( &d[12], vars->smp,  vars->swap );     // w4 SMP
    putFloat( &d[16], (float)hmt, vars->swap );     // w5 HMT
    putFloat( &d[20], (float)hmr, vars->swap );     // w6 HMR
    putInt(   &d[32], vars->nbyte, vars->swap );    // w9 NBYTE
    if( vars->nbyte == 1 ) {
      putFloat( &d[72], 0.0f, vars->swap );         // w19 ABIAS
      putInt(   &d[76], vars->lng2nrm, vars->swap );// w20 LNG2NRM
    }
    int w = 21;
    for( int ih = 0; ih < lhdr; ih++ ) {
      if( vars->hdrs[ih].spkName.empty() ) continue;
      char name[4] = { ' ', ' ', ' ', ' ' };
      memcpy( name, vars->hdrs[ih].spkName.c_str(), vars->hdrs[ih].spkName.size() );
      memcpy( &d[4*(w-1)], name, 4 );
      putFloat( &d[4*w], (float)(ih+1), vars->swap );
      w += 2;
    }
    memcpy( &d[4*(w-1)], "END ", 4 );
    putInt( &d[4*(nwords-1)], vars->ifw0, vars->swap );   // word IFW0-1 repeats IFW0 (as in ice-pak files)

    fseeko( vars->fp, 0, SEEK_SET );
    fwrite( &dir[0], 1, dir.size(), vars->fp );
    fclose( vars->fp );
    vars->fp = NULL;

    writer->line("");
    writer->line("  SEISPAK file written: %s", vars->filename.c_str());
    writer->line("     Traces: %ld  Traces/record (HMT): %d  Records (HMR): %d  Blocks: %lld", ntr, hmt, hmr, nblk);
    return true;
  }
}

using namespace mod_output_seispak;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_output_seispak_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->fp = NULL;
  vars->swap = false;
  vars->nbyte = 2;
  vars->lng2nrm = 100;
  vars->ntracesPerRecord = 0;
  vars->traceCounter = 0;
  vars->buffer = NULL;
  vars->isOK = false;

  std::string text;
  param->getString( "filename", &vars->filename );

  if( param->exists("format") ) {
    param->getString( "format", &text );
    if( !text.compare("int16") )      vars->nbyte = 2;
    else if( !text.compare("float") ) vars->nbyte = 4;
    else if( !text.compare("byte") )  vars->nbyte = 1;
    else writer->error("Unknown option for 'format': %s", text.c_str());
  }
  if( param->exists("lng2nrm") ) {
    param->getInt( "lng2nrm", &vars->lng2nrm );
    if( vars->lng2nrm <= 0 ) writer->error("Parameter 'lng2nrm' must be > 0");
  }
  bool bigEndian = !isLittleEndianMachine();
  if( param->exists("byte_order") ) {
    param->getString( "byte_order", &text );
    if( !text.compare("little") )   bigEndian = false;
    else if( !text.compare("big") ) bigEndian = true;
    else writer->error("Unknown option for 'byte_order': %s", text.c_str());
  }
  vars->swap = ( bigEndian == isLittleEndianMachine() );

  if( param->exists("ntraces_record") ) {
    param->getInt( "ntraces_record", &vars->ntracesPerRecord );
    if( vars->ntracesPerRecord < 0 ) vars->ntracesPerRecord = 0;
  }

  //---------------------------------------------
  // SEISPAK trace headers
  int nLines = param->getNumLines( "header" );
  if( nLines > 0 ) {
    for( int il = 0; il < nLines; il++ ) {
      HdrOut h;
      std::string cname;
      param->getStringAtLine( "header", &cname, il, 0 );
      std::string sname = "";
      if( param->getNumValues( "header", il ) > 1 ) param->getStringAtLine( "header", &sname, il, 1 );
      if( !cname.compare("_record") ) { h.source = SRC_RECORD; h.hdrId = -1; }
      else if( !cname.compare("_trace") ) { h.source = SRC_TRACE; h.hdrId = -1; }
      else {
        if( !hdef->headerExists( cname ) ) writer->error("Trace header '%s' does not exist", cname.c_str());
        type_t t = hdef->headerType( cname );
        if( t == TYPE_STRING ) writer->error("Trace header '%s' is a string header; numeric header required", cname.c_str());
        h.source = SRC_HEADER;
        h.hdrId = hdef->headerIndex( cname );
      }
      if( sname.empty() ) {
        std::string base = cname;
        if( base.compare(0,4,"spk_") == 0 ) base = base.substr(4);
        if( base[0] == '_' ) base = base.substr(1);
        sname = base;
      }
      if( !sname.compare("-") ) sname = "";
      h.spkName = toUpper4( sname );
      vars->hdrs.push_back( h );
    }
  }
  else {
    // Automatic: re-use SEISPAK headers created by INPUT_SEISPAK (in order of definition)
    for( int i = 0; i < hdef->numHeaders(); i++ ) {
      std::string name = hdef->headerName( i );
      if( name.compare(0,4,"spk_") != 0 ) continue;
      if( !name.compare("spk_rec") || !name.compare("spk_trc") ) continue;
      HdrOut h;
      h.source = SRC_HEADER;
      h.hdrId  = i;
      std::string base = name.substr(4);
      bool generic = ( base.size() > 1 && base[0] == 'h' && base.find_first_not_of("0123456789",1) == std::string::npos );
      h.spkName = generic ? "" : toUpper4( base );
      vars->hdrs.push_back( h );
    }
    if( vars->hdrs.empty() ) {
      HdrOut r; r.source = SRC_RECORD; r.hdrId = -1; r.spkName = "RNUM";
      HdrOut t; t.source = SRC_TRACE;  t.hdrId = -1; t.spkName = "TNUM";
      vars->hdrs.push_back( r );
      vars->hdrs.push_back( t );
    }
  }
  int lhdr = (int)vars->hdrs.size();
  int nNamed = 0;
  for( int ih = 0; ih < lhdr; ih++ ) if( !vars->hdrs[ih].spkName.empty() ) nNamed++;

  //---------------------------------------------
  vars->nt  = shdr->numSamples;
  vars->smp = shdr->sampleInt;
  vars->ifw0 = 21 + 2*nNamed + 2;     // words 21.. name/value pairs, 'END ', IFW0 repeat
  if( vars->nbyte == 4 ) {
    vars->recordBytes = 4LL * ( lhdr + vars->nt );
  }
  else if( vars->nbyte == 2 ) {
    vars->recordBytes = 4LL * ( lhdr + vars->nt/2 + ( (vars->nt % 2 == 0) ? 1 : 2 ) );
  }
  else {
    vars->recordBytes = 4LL*lhdr + vars->nt + 4LL*( (vars->nt + vars->lng2nrm - 1) / vars->lng2nrm );
  }
  vars->firstTraceByte = BASE_OFFSET + 4LL * ( vars->ifw0 - 1 );
  vars->buffer = new char[ (size_t)vars->recordBytes ];
  vars->work.resize( lhdr > 0 ? lhdr : 1 );

  vars->fp = fopen( vars->filename.c_str(), "w+b" );
  if( vars->fp == NULL ) writer->error("Cannot open output file '%s'", vars->filename.c_str());
  vars->isOK = true;

  char const* fmtName[5] = { "", "scaled byte", "scaled int16", "", "float32" };
  writer->line("");
  writer->line("  Output file:  %s", vars->filename.c_str());
  writer->line("  NBYTE: %d (%s), byte order: %s", vars->nbyte, fmtName[vars->nbyte], bigEndian ? "big endian" : "little endian");
  writer->line("  NT: %d  SMP: %f ms  LHDR: %d  IFW0: %d", vars->nt, vars->smp, lhdr, vars->ifw0);
  writer->line("  Traces/record: %s", vars->ntracesPerRecord > 0 ? "fixed" : "all traces in one record");
  for( int ih = 0; ih < lhdr; ih++ ) {
    std::string src = ( vars->hdrs[ih].source == SRC_RECORD ) ? "(record number)" :
                      ( vars->hdrs[ih].source == SRC_TRACE )  ? "(trace number in record)" :
                      hdef->headerName( vars->hdrs[ih].hdrId );
    writer->line("  Trace header word %2d  %-4s  <-  %s", ih+1,
                 vars->hdrs[ih].spkName.empty() ? "----" : vars->hdrs[ih].spkName.c_str(), src.c_str());
  }
  writer->line("");
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_output_seispak_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csTrace* trace = traceGather->trace(0);
  csTraceHeader* trcHdr = trace->getTraceHeader();
  float const* samples = trace->getTraceSamples();

  long k = vars->traceCounter;
  int hmt = vars->ntracesPerRecord;
  int lhdr = (int)vars->hdrs.size();
  for( int ih = 0; ih < lhdr; ih++ ) {
    HdrOut const& h = vars->hdrs[ih];
    if( h.source == SRC_HEADER )      vars->work[ih] = (float)trcHdr->doubleValue( h.hdrId );
    else if( h.source == SRC_RECORD ) vars->work[ih] = (float)( hmt > 0 ? k / hmt + 1 : 1 );
    else                              vars->work[ih] = (float)( hmt > 0 ? k % hmt + 1 : k + 1 );
  }
  encodeTrace( vars, &vars->work[0], samples, vars->buffer );
  off_t pos = (off_t)( vars->firstTraceByte + (long long)k * vars->recordBytes );
  if( fseeko( vars->fp, pos, SEEK_SET ) != 0 ||
      fwrite( vars->buffer, 1, (size_t)vars->recordBytes, vars->fp ) != (size_t)vars->recordBytes ) {
    writer->error("Error writing trace %ld to SEISPAK file '%s' (disk full?)", k+1, vars->filename.c_str());
  }
  vars->traceCounter += 1;
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_output_seispak_( csParamDef* pdef ) {
  pdef->setModule( "OUTPUT_SEISPAK", "Output data in SEISPAK format (Texaco ice-pak, 'ICEp' files)",
                   "Writes SEISPAK disk files as the ice-pak library does (routines PUT2/PUT1/FLT2FX). "
                   "By default, SEISPAK trace headers created by INPUT_SEISPAK (spk_*) are written back in their original order, "
                   "so that INPUT_SEISPAK -> OUTPUT_SEISPAK preserves the file structure. If no spk_* headers exist, "
                   "two headers are written: RNUM (record number) and TNUM (trace number within record). "
                   "The header directory is written when the flow finishes." );

  pdef->addParam( "filename", "Output file name", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "Output file name" );

  pdef->addParam( "format", "Sample format (header word 9, NBYTE)", NUM_VALUES_FIXED );
  pdef->addValue( "int16", VALTYPE_OPTION );
  pdef->addOption( "int16", "NBYTE=2: 16bit integer, scaled per trace (max amplitude = 32000)" );
  pdef->addOption( "float", "NBYTE=4: 32bit float (no loss of precision)" );
  pdef->addOption( "byte",  "NBYTE=1: 8bit integer, scaled per group of 'lng2nrm' samples" );

  pdef->addParam( "lng2nrm", "Samples per normalisation group (format 'byte' only)", NUM_VALUES_FIXED );
  pdef->addValue( "100", VALTYPE_NUMBER, "Number of samples" );

  pdef->addParam( "ntraces_record", "Number of traces per SEISPAK record (HMT)", NUM_VALUES_FIXED,
                  "Set 0 to write all traces into one single record. If the total number of traces is not a multiple of this value, the last record is filled with zero traces." );
  pdef->addValue( "0", VALTYPE_NUMBER, "Traces per record" );

  pdef->addParam( "byte_order", "Byte order of output file", NUM_VALUES_FIXED );
  pdef->addValue( "little", VALTYPE_OPTION );
  pdef->addOption( "little", "Little endian (Linux PC, DEC)" );
  pdef->addOption( "big", "Big endian (SGI, SUN, HP, IBM)" );

  pdef->addParam( "header", "SEISPAK trace header word (one line per word, in order)", NUM_VALUES_VARIABLE,
                  "Overrides the automatic header selection. Use '_record' / '_trace' for the SEISPAK record / trace number." );
  pdef->addValue( "", VALTYPE_STRING, "SeaSeis trace header name, or _record / _trace" );
  pdef->addValue( "", VALTYPE_STRING, "SEISPAK header name (up to 4 characters) stored in the file's global list. Use '-' for an unnamed word." );
}

//************************************************************************************************
// Start exec phase
//*************************************************************************************************
bool start_exec_mod_output_seispak_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

//************************************************************************************************
// Cleanup phase
//*************************************************************************************************
void cleanup_mod_output_seispak_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( vars->fp != NULL ) {
    if( vars->isOK && vars->traceCounter > 0 ) {
      finalizeFile( vars, writer );
    }
    else {
      fclose( vars->fp );
      vars->fp = NULL;
    }
  }
  if( vars->buffer != NULL ) { delete [] vars->buffer; vars->buffer = NULL; }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_output_seispak_( csParamDef* pdef ) {
  params_mod_output_seispak_( pdef );
}
extern "C" void _init_mod_output_seispak_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_output_seispak_( param, env, writer );
}
extern "C" bool _start_exec_mod_output_seispak_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_output_seispak_( env, writer );
}
extern "C" void _exec_mod_output_seispak_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_output_seispak_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_output_seispak_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_output_seispak_( env, writer );
}
