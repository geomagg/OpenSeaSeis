/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/*
 * Module OUTPUT_HDF5: write traces to an HDF5 file in the OBN layout read by INPUT_HDF5.
 * Author: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026.
 *
 * Output layout:
 *   /data                   float  [nReceivers][nChannels][nSamples]  (chunked, one chunk per trace)
 *   /header/dt              double sample interval [s]
 *   /header/time            double start time, seconds since 1970-01-01 (UTC)   (if time headers exist)
 *   /header/channel_number  int    original value of the channel header for each channel index
 *   /metadata/REC_X         double receiver X
 *   /metadata/REC_Y         double receiver Y
 *   /metadata/elevation     double receiver elevation (negative below sea level)
 *   /metadata/sensor_id     string (variable length)
 *
 * Receivers are identified by the value(s) of the header(s) given in 'receiver_hdr' and channels by
 * the value of 'channel_hdr'. Both are numbered in order of first appearance in the flow.
 * Receiver/channel cells that receive no trace are filled with zeros.
 */

#include "cseis_includes.h"
#include <hdf5.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <vector>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

namespace mod_output_hdf5 {
  struct VariableStruct {
    std::string filename;
    int   compression;
    int   nSamp;
    double dt_s;

    hid_t fileId;
    hid_t dataId;
    hsize_t nRecAlloc;
    hsize_t nChanAlloc;

    std::vector<int> hdrId_recKey;
    int hdrId_chan;
    int hdrId_rec_x, hdrId_rec_y, hdrId_rec_elev, hdrId_rec_z, hdrId_node, hdrId_rec_line, hdrId_rcv, hdrId_sensor_id;
    int hdrId_time_samp1, hdrId_time_samp1_us;

    std::map<std::string,int> recMap;
    std::map<int,int> chanMap;
    std::vector<int> chanNumber;
    std::vector<double> recX, recY, elev;
    std::vector<std::string> sensorId;
    std::vector< std::vector<char> > written;   // [rec][chan]

    double time_s;          // start time of first trace (<0: unknown)
    int    numTimeDiffer;   // traces whose start time differs from the first one
    int    traceCounter;
  };

  static int hdrIdOrMinus1( csTraceHeaderDef const* hdef, char const* name ) {
    return hdef->headerExists( name ) ? hdef->headerIndex( name ) : -1;
  }

  static bool writeDoubles( hid_t loc, char const* name, std::vector<double> const& v ) {
    if( v.empty() ) return true;
    hsize_t n = v.size();
    hid_t sp = H5Screate_simple( 1, &n, NULL );
    hid_t ds = H5Dcreate2( loc, name, H5T_IEEE_F64LE, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    herr_t st = ( ds >= 0 ) ? H5Dwrite( ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0] ) : -1;
    if( ds >= 0 ) H5Dclose( ds );
    H5Sclose( sp );
    return st >= 0;
  }

  static bool writeInts( hid_t loc, char const* name, std::vector<int> const& v ) {
    if( v.empty() ) return true;
    hsize_t n = v.size();
    hid_t sp = H5Screate_simple( 1, &n, NULL );
    hid_t ds = H5Dcreate2( loc, name, H5T_STD_I32LE, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    herr_t st = ( ds >= 0 ) ? H5Dwrite( ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v[0] ) : -1;
    if( ds >= 0 ) H5Dclose( ds );
    H5Sclose( sp );
    return st >= 0;
  }

  static bool writeStrings( hid_t loc, char const* name, std::vector<std::string> const& v ) {
    if( v.empty() ) return true;
    hsize_t n = v.size();
    std::vector<char const*> ptr( v.size() );
    for( size_t i = 0; i < v.size(); i++ ) ptr[i] = v[i].c_str();
    hid_t type = H5Tcopy( H5T_C_S1 );
    H5Tset_size( type, H5T_VARIABLE );
    H5Tset_cset( type, H5T_CSET_UTF8 );
    hid_t sp = H5Screate_simple( 1, &n, NULL );
    hid_t ds = H5Dcreate2( loc, name, type, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    herr_t st = ( ds >= 0 ) ? H5Dwrite( ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ptr[0] ) : -1;
    if( ds >= 0 ) H5Dclose( ds );
    H5Sclose( sp );
    H5Tclose( type );
    return st >= 0;
  }

  // Close and delete an incomplete output file (called before a fatal error in the exec phase)
  static void abortFile( VariableStruct* vars ) {
    if( vars->dataId >= 0 ) { H5Dclose( vars->dataId ); vars->dataId = -1; }
    if( vars->fileId >= 0 ) { H5Fclose( vars->fileId ); vars->fileId = -1; }
    remove( vars->filename.c_str() );
  }

  // Write header + metadata groups and close the file
  static void finishFile( VariableStruct* vars, csLogWriter* writer ) {
    if( vars->fileId < 0 ) return;
    if( vars->dataId >= 0 ) { H5Dclose( vars->dataId ); vars->dataId = -1; }

    bool ok = true;
    hid_t g = H5Gcreate2( vars->fileId, "/header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    ok &= writeDoubles( g, "dt", std::vector<double>( 1, vars->dt_s ) );
    if( vars->time_s >= 0 ) ok &= writeDoubles( g, "time", std::vector<double>( 1, vars->time_s ) );
    ok &= writeInts( g, "channel_number", vars->chanNumber );
    H5Gclose( g );

    g = H5Gcreate2( vars->fileId, "/metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
    if( vars->hdrId_rec_x >= 0 ) ok &= writeDoubles( g, "REC_X", vars->recX );
    if( vars->hdrId_rec_y >= 0 ) ok &= writeDoubles( g, "REC_Y", vars->recY );
    if( vars->hdrId_rec_elev >= 0 || vars->hdrId_rec_z >= 0 ) ok &= writeDoubles( g, "elevation", vars->elev );
    bool anyId = false;
    for( size_t i = 0; i < vars->sensorId.size(); i++ ) if( !vars->sensorId[i].empty() ) { anyId = true; break; }
    if( anyId ) ok &= writeStrings( g, "sensor_id", vars->sensorId );
    H5Gclose( g );

    H5Fclose( vars->fileId );
    vars->fileId = -1;

    if( writer != NULL ) {
      int nEmpty = 0;
      for( size_t r = 0; r < vars->written.size(); r++ )
        for( size_t c = 0; c < vars->chanNumber.size(); c++ )
          if( c >= vars->written[r].size() || !vars->written[r][c] ) nEmpty++;
      writer->line( "OUTPUT_HDF5: file '%s' closed", vars->filename.c_str() );
      writer->line( "  traces written: %d   receivers: %d   channels: %d   samples: %d   dt: %g s",
                    vars->traceCounter, (int)vars->recMap.size(), (int)vars->chanNumber.size(), vars->nSamp, vars->dt_s );
      if( nEmpty > 0 ) writer->line( "  NOTE: %d receiver/channel cells received no trace and are filled with zeros", nEmpty );
      if( vars->numTimeDiffer > 0 ) writer->line( "  WARNING: %d traces have a start time different from the first trace (only one start time is stored in /header/time)", vars->numTimeDiffer );
      if( !ok ) writer->line( "  WARNING: error while writing header/metadata data sets" );
    }
  }
}
using namespace mod_output_hdf5;

//*************************************************************************************************
// Init phase
//
void init_mod_output_hdf5_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->fileId = vars->dataId = -1;
  vars->nRecAlloc = vars->nChanAlloc = 0;
  vars->compression = 0;
  vars->time_s = -1.0;
  vars->numTimeDiffer = 0;
  vars->traceCounter = 0;
  vars->nSamp = shdr->numSamples;
  vars->dt_s  = (double)shdr->sampleInt / 1000.0;

  H5Eset_auto2( H5E_DEFAULT, NULL, NULL );

  //--- Parameters
  param->getString( "filename", &vars->filename );

  std::vector<std::string> recHdrNames;
  bool recAuto = !param->exists("receiver_hdr");
  if( !recAuto ) {
    int n = param->getNumValues( "receiver_hdr" );
    for( int i = 0; i < n; i++ ) {
      std::string name;
      param->getString( "receiver_hdr", &name, i );
      recHdrNames.push_back( toLowerCase( name ) );
    }
  }
  else {
    // Automatic: first available of rec_index / node / rec_line+rcv / rec_x+rec_y
    char const* cand[4][2] = { {"rec_index",NULL}, {"node",NULL}, {"rec_line","rcv"}, {"rec_x","rec_y"} };
    for( int k = 0; k < 4 && recHdrNames.empty(); k++ ) {
      if( hdef->headerExists( cand[k][0] ) && ( cand[k][1] == NULL || hdef->headerExists( cand[k][1] ) ) ) {
        recHdrNames.push_back( cand[k][0] );
        if( cand[k][1] != NULL ) recHdrNames.push_back( cand[k][1] );
      }
    }
    if( recHdrNames.empty() ) writer->error("No receiver header found (tried rec_index, node, rec_line+rcv, rec_x+rec_y). Specify parameter 'receiver_hdr'.");
  }
  for( size_t i = 0; i < recHdrNames.size(); i++ ) {
    if( !hdef->headerExists( recHdrNames[i] ) ) {
      writer->error("Receiver header '%s' does not exist. Set parameter 'receiver_hdr' to the header(s) identifying the receiver, e.g. 'node', 'rec_line rcv' or 'rec_x rec_y'.", recHdrNames[i].c_str());
    }
    if( hdef->headerType( recHdrNames[i] ) == TYPE_STRING ) writer->error("Receiver header '%s' is a string header. Use numeric headers.", recHdrNames[i].c_str());
    vars->hdrId_recKey.push_back( hdef->headerIndex( recHdrNames[i] ) );
  }

  std::string chanHdrName = "chan";
  if( param->exists("channel_hdr") ) {
    param->getString( "channel_hdr", &chanHdrName );
    chanHdrName = toLowerCase( chanHdrName );
  }
  if( !hdef->headerExists( chanHdrName ) ) writer->error("Channel header '%s' does not exist (parameter 'channel_hdr').", chanHdrName.c_str());
  vars->hdrId_chan = hdef->headerIndex( chanHdrName );

  if( param->exists("compression") ) {
    param->getInt( "compression", &vars->compression );
    if( vars->compression < 0 || vars->compression > 9 ) writer->error("Compression level must be between 0 (none) and 9.");
    if( vars->compression > 0 && H5Zfilter_avail( H5Z_FILTER_DEFLATE ) <= 0 ) {
      writer->warning("gzip compression not available in this HDF5 library. Writing uncompressed.");
      vars->compression = 0;
    }
  }

  //--- Optional headers used for metadata
  vars->hdrId_rec_x    = hdrIdOrMinus1( hdef, "rec_x" );
  vars->hdrId_rec_y    = hdrIdOrMinus1( hdef, "rec_y" );
  vars->hdrId_rec_elev = hdrIdOrMinus1( hdef, "rec_elev" );
  vars->hdrId_rec_z    = hdrIdOrMinus1( hdef, "rec_z" );
  vars->hdrId_node     = hdrIdOrMinus1( hdef, "node" );
  vars->hdrId_rec_line = hdrIdOrMinus1( hdef, "rec_line" );
  vars->hdrId_rcv      = hdrIdOrMinus1( hdef, "rcv" );
  vars->hdrId_sensor_id = -1;
  if( hdef->headerExists( "sensor_id" ) && hdef->headerType( "sensor_id" ) == TYPE_STRING ) vars->hdrId_sensor_id = hdef->headerIndex( "sensor_id" );
  vars->hdrId_time_samp1    = hdrIdOrMinus1( hdef, "time_samp1" );
  vars->hdrId_time_samp1_us = hdrIdOrMinus1( hdef, "time_samp1_us" );

  //--- Create file and extendible data set [rec][chan][sample]
  vars->fileId = H5Fcreate( vars->filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
  if( vars->fileId < 0 ) writer->error("Cannot create HDF5 file '%s'.", vars->filename.c_str());

  hsize_t dims[3]    = { 0, 0, (hsize_t)vars->nSamp };
  hsize_t maxdims[3] = { H5S_UNLIMITED, H5S_UNLIMITED, (hsize_t)vars->nSamp };
  hsize_t chunk[3]   = { 1, 1, (hsize_t)vars->nSamp };
  hid_t space = H5Screate_simple( 3, dims, maxdims );
  hid_t plist = H5Pcreate( H5P_DATASET_CREATE );
  H5Pset_chunk( plist, 3, chunk );
  float fill = 0.0f;
  H5Pset_fill_value( plist, H5T_NATIVE_FLOAT, &fill );
  if( vars->compression > 0 ) {
    H5Pset_shuffle( plist );
    H5Pset_deflate( plist, (unsigned)vars->compression );
  }
  vars->dataId = H5Dcreate2( vars->fileId, "/data", H5T_IEEE_F32LE, space, H5P_DEFAULT, plist, H5P_DEFAULT );
  H5Pclose( plist );
  H5Sclose( space );
  if( vars->dataId < 0 ) writer->error("Cannot create data set '/data' in file '%s'.", vars->filename.c_str());

  //--- Log
  std::string keyText;
  for( size_t i = 0; i < recHdrNames.size(); i++ ) keyText += ( i ? " + " : "" ) + recHdrNames[i];
  writer->line( "  Output file:      %s", vars->filename.c_str() );
  writer->line( "  Receiver key:     %s%s", keyText.c_str(), recAuto ? "  (automatic)" : "" );
  writer->line( "  Channel key:      %s", chanHdrName.c_str() );
  writer->line( "  Samples / dt:     %d / %g s", vars->nSamp, vars->dt_s );
  writer->line( "  Compression:      %s", vars->compression > 0 ? "gzip + shuffle" : "none" );
  if( vars->hdrId_rec_x < 0 || vars->hdrId_rec_y < 0 ) writer->line( "  NOTE: header rec_x/rec_y not found: /metadata/REC_X, REC_Y not written" );
  if( vars->hdrId_time_samp1 < 0 ) writer->line( "  NOTE: header time_samp1 not found: /header/time not written" );
}

//*************************************************************************************************
// Exec phase
//
void exec_mod_output_hdf5_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );

  if( vars->dataId < 0 ) return;
  csTrace* trace = traceGather->trace(0);
  csTraceHeader const* trcHdr = trace->getTraceHeader();

  //--- Receiver index (order of first appearance)
  std::string key;
  char buf[64];
  for( size_t i = 0; i < vars->hdrId_recKey.size(); i++ ) {
    snprintf( buf, sizeof(buf), "%.15g|", trcHdr->doubleValue( vars->hdrId_recKey[i] ) );
    key += buf;
  }
  int irec;
  std::map<std::string,int>::iterator itRec = vars->recMap.find( key );
  if( itRec == vars->recMap.end() ) {
    irec = (int)vars->recMap.size();
    vars->recMap[key] = irec;
    // Per receiver metadata, taken from the first trace of this receiver
    vars->recX.push_back( vars->hdrId_rec_x >= 0 ? trcHdr->doubleValue( vars->hdrId_rec_x ) : 0.0 );
    vars->recY.push_back( vars->hdrId_rec_y >= 0 ? trcHdr->doubleValue( vars->hdrId_rec_y ) : 0.0 );
    double e = 0.0;
    if( vars->hdrId_rec_elev >= 0 ) e = trcHdr->doubleValue( vars->hdrId_rec_elev );
    else if( vars->hdrId_rec_z >= 0 ) e = trcHdr->doubleValue( vars->hdrId_rec_z );
    vars->elev.push_back( e );
    std::string id;
    if( vars->hdrId_sensor_id >= 0 ) id = trcHdr->stringValue( vars->hdrId_sensor_id );
    size_t z = id.find( '\0' );
    if( z != std::string::npos ) id = id.substr( 0, z );
    while( !id.empty() && id[id.size()-1] == ' ' ) id.erase( id.size()-1 );
    if( id.empty() && vars->hdrId_node >= 0 && trcHdr->intValue( vars->hdrId_node ) > 0 ) {
      snprintf( buf, sizeof(buf), "RID%d", trcHdr->intValue( vars->hdrId_node ) );
      id = buf;
    }
    else if( id.empty() && vars->hdrId_rec_line >= 0 && vars->hdrId_rcv >= 0 &&
             ( trcHdr->intValue( vars->hdrId_rec_line ) > 0 || trcHdr->intValue( vars->hdrId_rcv ) > 0 ) ) {
      snprintf( buf, sizeof(buf), "RID%04d%04d", trcHdr->intValue( vars->hdrId_rec_line ), trcHdr->intValue( vars->hdrId_rcv ) );
      id = buf;
    }
    vars->sensorId.push_back( id );
    vars->written.push_back( std::vector<char>() );
  }
  else {
    irec = itRec->second;
  }

  //--- Channel index (order of first appearance)
  int chanValue = trcHdr->intValue( vars->hdrId_chan );
  int ichan;
  std::map<int,int>::iterator itChan = vars->chanMap.find( chanValue );
  if( itChan == vars->chanMap.end() ) {
    ichan = (int)vars->chanMap.size();
    vars->chanMap[chanValue] = ichan;
    vars->chanNumber.push_back( chanValue );
  }
  else {
    ichan = itChan->second;
  }

  std::vector<char>& wr = vars->written[irec];
  if( (int)wr.size() <= ichan ) wr.resize( ichan+1, 0 );
  if( wr[ichan] ) {
    abortFile( vars );
    writer->error("Trace %d: receiver #%d / channel %d was already written. Two input traces map to the same receiver and channel. "
                  "Use a receiver key that is unique per receiver (parameter 'receiver_hdr', e.g. 'fileno rec_index' when several input files are combined).",
                  vars->traceCounter+1, irec+1, chanValue );
  }
  wr[ichan] = 1;

  //--- Start time
  if( vars->hdrId_time_samp1 >= 0 ) {
    double t = (double)trcHdr->intValue( vars->hdrId_time_samp1 );
    if( vars->hdrId_time_samp1_us >= 0 ) t += 1.0e-6 * (double)trcHdr->intValue( vars->hdrId_time_samp1_us );
    if( vars->traceCounter == 0 ) vars->time_s = t;
    else if( fabs( t - vars->time_s ) > 1.0e-6 ) vars->numTimeDiffer++;
  }

  //--- Extend data set if needed, then write one trace
  hsize_t needRec  = std::max( vars->nRecAlloc,  (hsize_t)irec+1 );
  hsize_t needChan = std::max( vars->nChanAlloc, (hsize_t)ichan+1 );
  if( needRec != vars->nRecAlloc || needChan != vars->nChanAlloc ) {
    if( needRec > vars->nRecAlloc ) needRec = std::max( needRec, vars->nRecAlloc + vars->nRecAlloc/2 + 16 );  // grow in blocks
    hsize_t newDims[3] = { needRec, needChan, (hsize_t)vars->nSamp };
    if( H5Dset_extent( vars->dataId, newDims ) < 0 ) { abortFile( vars ); writer->error("Cannot extend data set '/data'."); }
    vars->nRecAlloc = needRec;
    vars->nChanAlloc = needChan;
  }
  hid_t fspace = H5Dget_space( vars->dataId );
  hsize_t start[3] = { (hsize_t)irec, (hsize_t)ichan, 0 };
  hsize_t count[3] = { 1, 1, (hsize_t)vars->nSamp };
  H5Sselect_hyperslab( fspace, H5S_SELECT_SET, start, NULL, count, NULL );
  hid_t mspace = H5Screate_simple( 1, &count[2], NULL );
  herr_t st = H5Dwrite( vars->dataId, H5T_NATIVE_FLOAT, mspace, fspace, H5P_DEFAULT, trace->getTraceSamples() );
  H5Sclose( mspace );
  H5Sclose( fspace );
  if( st < 0 ) {
    int n = vars->traceCounter+1;
    abortFile( vars );
    writer->error("Error writing trace %d to file '%s' (disk full?). Incomplete file deleted.", n, vars->filename.c_str());
  }

  vars->traceCounter++;
}

//*************************************************************************************************
// Parameter definition (help: seaseis -m output_hdf5)
//
void params_mod_output_hdf5_( csParamDef* pdef ) {
  pdef->setModule( "OUTPUT_HDF5", "Output OBN data to HDF5 file",
                   "Writes traces to data set /data [receiver][channel][sample], plus /header and /metadata, in the layout read by INPUT_HDF5" );

  pdef->addDoc("Receivers are identified by the header(s) in 'receiver_hdr', channels by 'channel_hdr'; both are numbered in order of first appearance.");
  pdef->addDoc("Per receiver: /metadata/REC_X, REC_Y (rec_x, rec_y), elevation (rec_elev, or rec_z), sensor_id (header sensor_id, else 'RID'+node, else 'RID'+rec_line+rcv).");
  pdef->addDoc("/header/dt [s], /header/time (time_samp1 + time_samp1_us of the first trace, UTC) and /header/channel_number (value of channel_hdr per channel).");
  pdef->addDoc("Receiver/channel cells without a trace are filled with zeros. Two traces with the same receiver and channel are an error (the incomplete file is deleted).");

  pdef->addParam( "filename", "Output HDF5 file name", NUM_VALUES_FIXED, "An existing file is overwritten" );
  pdef->addValue( "", VALTYPE_STRING, "Output file name" );

  pdef->addParam( "receiver_hdr", "Trace header(s) identifying the receiver", NUM_VALUES_VARIABLE,
                  "Default: first available of rec_index, node, rec_line+rcv, rec_x+rec_y. Several headers are combined, e.g. 'fileno rec_index'" );
  pdef->addValue( "", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "channel_hdr", "Trace header identifying the channel", NUM_VALUES_FIXED );
  pdef->addValue( "chan", VALTYPE_STRING, "Trace header name" );

  pdef->addParam( "compression", "gzip compression level", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "0: no compression, 1 (fast) - 9 (smallest file)" );
}

//************************************************************************************************
// Start exec phase
//
bool start_exec_mod_output_hdf5_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

//************************************************************************************************
// Cleanup phase
//
void cleanup_mod_output_hdf5_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( vars->dataId >= 0 && vars->nRecAlloc != vars->recMap.size() ) {
    // Trim receiver dimension to the number of receivers actually written
    hsize_t dims[3] = { (hsize_t)vars->recMap.size(), vars->nChanAlloc, (hsize_t)vars->nSamp };
    H5Dset_extent( vars->dataId, dims );
  }
  finishFile( vars, writer );
  delete vars; vars = NULL;
}

extern "C" void _params_mod_output_hdf5_( csParamDef* pdef ) {
  params_mod_output_hdf5_( pdef );
}
extern "C" void _init_mod_output_hdf5_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_output_hdf5_( param, env, writer );
}
extern "C" bool _start_exec_mod_output_hdf5_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_output_hdf5_( env, writer );
}
extern "C" void _exec_mod_output_hdf5_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_output_hdf5_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_output_hdf5_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_output_hdf5_( env, writer );
}
