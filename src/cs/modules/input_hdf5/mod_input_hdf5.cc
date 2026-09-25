/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/*
 * Module INPUT_HDF5: read OBN data stored in HDF5 files.
 * Author: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026.
 *
 * Expected layout (dataset names can be changed by user parameters):
 *   /data                float  [nReceivers][nChannels][nSamples]
 *   /header/dt           double sample interval [s]
 *   /header/time         double start time of the record, seconds since 1970-01-01 (UTC)
 *   /metadata/REC_X      receiver X          (optional)
 *   /metadata/REC_Y      receiver Y          (optional)
 *   /metadata/elevation  receiver elevation  (optional, negative below sea level)
 *   /metadata/sensor_id  string, e.g. "..._RID116117011.P26774" (optional)
 *                        digits after 'RID': 1-4 = receiver line, 5-8 = receiver station
 *
 * Each (receiver, channel) pair becomes one output trace.
 */

#include "cseis_includes.h"
#include <hdf5.h>
#include <cmath>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

namespace mod_input_hdf5 {
  static int const ORDER_CHANNEL  = 1;
  static int const ORDER_RECEIVER = 2;

  struct FileInfo {
    std::string name;
    hsize_t nRec;
    hsize_t nChan;
    hsize_t nSamp;
    double  dt_s;
    double  time_s;          // start time, seconds since 1970 (UTC); < 0 if not available
    std::vector<double> recX, recY, elev;
    std::vector<std::string> sensorId;
    std::vector<double> chanNumber;   // optional /header/channel_number (written by OUTPUT_HDF5)
  };

  struct VariableStruct {
    std::vector<FileInfo> files;
    std::vector<int> channels;     // 0-based channel indices to output
    int order;
    std::string dsData;
    int ridLine1, ridLine2, ridSta1, ridSta2;   // 1-based digit positions in RID number

    // Reading state
    int   fileIndex;
    hid_t fileId;
    hid_t dataId;
    hid_t spaceId;
    hsize_t recIndex;
    int     chanPos;               // position in 'channels'
    int     traceCounter;
    bool    atEOF;

    // Header indices
    int hdrId_trcno, hdrId_fileno, hdrId_chan, hdrId_node, hdrId_rec_line, hdrId_rcv, hdrId_rec_index;
    int hdrId_rec_x, hdrId_rec_y, hdrId_rec_z, hdrId_rec_elev, hdrId_rec_wdep;
    int hdrId_time_samp1, hdrId_time_samp1_us, hdrId_time_code;
    int hdrId_time_year, hdrId_time_day, hdrId_time_hour, hdrId_time_min, hdrId_time_sec;
    int hdrId_sensor_id;
  };

  //--------------------------------------------------------------------
  static bool linkExists( hid_t fid, std::string const& path ) {
    // Check each component of the path (H5Lexists needs parent groups to exist)
    size_t pos = 1;
    while( pos <= path.size() ) {
      size_t next = path.find( '/', pos );
      std::string part = path.substr( 0, next == std::string::npos ? path.size() : next );
      if( H5Lexists( fid, part.c_str(), H5P_DEFAULT ) <= 0 ) return false;
      if( next == std::string::npos ) break;
      pos = next + 1;
    }
    return true;
  }

  static bool readNumbers( hid_t fid, std::string const& path, std::vector<double>& out ) {
    out.clear();
    if( !linkExists( fid, path ) ) return false;
    hid_t ds = H5Dopen2( fid, path.c_str(), H5P_DEFAULT );
    if( ds < 0 ) return false;
    hid_t sp = H5Dget_space( ds );
    hssize_t n = H5Sget_simple_extent_npoints( sp );
    bool ok = false;
    if( n > 0 ) {
      out.resize( (size_t)n );
      ok = H5Dread( ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &out[0] ) >= 0;
    }
    H5Sclose( sp );
    H5Dclose( ds );
    if( !ok ) out.clear();
    return ok;
  }

  static bool readStrings( hid_t fid, std::string const& path, std::vector<std::string>& out ) {
    out.clear();
    if( !linkExists( fid, path ) ) return false;
    hid_t ds = H5Dopen2( fid, path.c_str(), H5P_DEFAULT );
    if( ds < 0 ) return false;
    hid_t ftype = H5Dget_type( ds );
    hid_t sp = H5Dget_space( ds );
    hssize_t n = H5Sget_simple_extent_npoints( sp );
    bool ok = false;
    if( n > 0 && H5Tget_class( ftype ) == H5T_STRING ) {
      if( H5Tis_variable_str( ftype ) > 0 ) {
        std::vector<char*> buf( (size_t)n, (char*)NULL );
        hid_t mtype = H5Tcopy( H5T_C_S1 );
        H5Tset_size( mtype, H5T_VARIABLE );
        H5Tset_cset( mtype, H5Tget_cset( ftype ) );
        if( H5Dread( ds, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buf[0] ) >= 0 ) {
          for( size_t i = 0; i < buf.size(); i++ ) out.push_back( buf[i] ? std::string( buf[i] ) : std::string() );
          H5Dvlen_reclaim( mtype, sp, H5P_DEFAULT, &buf[0] );
          ok = true;
        }
        H5Tclose( mtype );
      }
      else {
        size_t len = H5Tget_size( ftype );
        std::vector<char> buf( (size_t)n * len + 1, 0 );
        hid_t mtype = H5Tcopy( H5T_C_S1 );
        H5Tset_size( mtype, len );
        if( H5Dread( ds, mtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buf[0] ) >= 0 ) {
          for( hssize_t i = 0; i < n; i++ ) {
            std::string s( &buf[i*len], len );
            size_t z = s.find( '\0' );
            out.push_back( z == std::string::npos ? s : s.substr( 0, z ) );
          }
          ok = true;
        }
        H5Tclose( mtype );
      }
    }
    H5Sclose( sp );
    H5Tclose( ftype );
    H5Dclose( ds );
    if( !ok ) out.clear();
    return ok;
  }

  // Digits following 'RID' in the sensor id, or empty string
  static std::string ridDigits( std::string const& id ) {
    size_t p = id.find( "RID" );
    if( p == std::string::npos ) return std::string();
    p += 3;
    size_t e = p;
    while( e < id.size() && id[e] >= '0' && id[e] <= '9' ) e++;
    return id.substr( p, e - p );
  }

  static int digitsToInt( std::string const& digits, int pos1, int pos2 ) {
    if( pos1 < 1 || pos2 < pos1 || (int)digits.size() < pos2 ) return 0;
    return atoi( digits.substr( pos1-1, pos2-pos1+1 ).c_str() );
  }
}
using namespace mod_input_hdf5;

//*************************************************************************************************
// Init phase
//
void init_mod_input_hdf5_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csTraceHeaderDef* hdef = env->headerDef;
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );

  edef->setExecType( EXEC_TYPE_INPUT );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->order    = ORDER_CHANNEL;
  vars->dsData   = "/data";
  vars->ridLine1 = 1; vars->ridLine2 = 4;
  vars->ridSta1  = 5; vars->ridSta2  = 8;
  vars->fileIndex = -1;
  vars->fileId = vars->dataId = vars->spaceId = -1;
  vars->recIndex = 0;
  vars->chanPos  = 0;
  vars->traceCounter = 0;
  vars->atEOF = false;

  H5Eset_auto2( H5E_DEFAULT, NULL, NULL );   // Errors are reported by this module

  //--- Parameters
  std::string dsDt = "/header/dt", dsTime = "/header/time", grpMeta = "/metadata";
  if( param->exists("dataset") )   param->getString( "dataset", &vars->dsData );
  if( param->exists("ds_dt") )     param->getString( "ds_dt", &dsDt );
  if( param->exists("ds_time") )   param->getString( "ds_time", &dsTime );
  if( param->exists("metadata") )  param->getString( "metadata", &grpMeta );

  if( param->exists("order") ) {
    std::string text;
    param->getString( "order", &text );
    text = toLowerCase( text );
    if( !text.compare("channel") )       vars->order = ORDER_CHANNEL;
    else if( !text.compare("receiver") ) vars->order = ORDER_RECEIVER;
    else writer->error("Unknown option for parameter 'order': '%s'", text.c_str());
  }
  if( param->exists("rid_line") ) {
    param->getInt( "rid_line", &vars->ridLine1, 0 );
    param->getInt( "rid_line", &vars->ridLine2, 1 );
  }
  if( param->exists("rid_station") ) {
    param->getInt( "rid_station", &vars->ridSta1, 0 );
    param->getInt( "rid_station", &vars->ridSta2, 1 );
  }

  //--- Input files (one or more 'filename' lines, each with one or more names)
  int nLines = param->getNumLines( "filename" );
  for( int iline = 0; iline < nLines; iline++ ) {
    int nValues = param->getNumValues( "filename", iline );
    for( int ival = 0; ival < nValues; ival++ ) {
      std::string name;
      param->getStringAtLine( "filename", &name, iline, ival );
      FileInfo info;
      info.name = name;
      vars->files.push_back( info );
    }
  }
  if( vars->files.empty() ) writer->error("No input file specified (parameter 'filename').");

  //--- Scan all files: dimensions, sample interval, metadata
  for( size_t ifile = 0; ifile < vars->files.size(); ifile++ ) {
    FileInfo& f = vars->files[ifile];
    hid_t fid = H5Fopen( f.name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
    if( fid < 0 ) writer->error("Cannot open HDF5 file '%s'.", f.name.c_str());
    if( !linkExists( fid, vars->dsData ) ) writer->error("Dataset '%s' not found in file '%s'.", vars->dsData.c_str(), f.name.c_str());
    hid_t ds = H5Dopen2( fid, vars->dsData.c_str(), H5P_DEFAULT );
    hid_t sp = H5Dget_space( ds );
    int ndims = H5Sget_simple_extent_ndims( sp );
    if( ndims != 3 ) writer->error("Dataset '%s' in file '%s' has %d dimensions. Expected 3: [receiver][channel][sample].", vars->dsData.c_str(), f.name.c_str(), ndims);
    hsize_t dims[3];
    H5Sget_simple_extent_dims( sp, dims, NULL );
    f.nRec = dims[0]; f.nChan = dims[1]; f.nSamp = dims[2];
    H5Sclose( sp );
    H5Dclose( ds );

    std::vector<double> v;
    if( !readNumbers( fid, dsDt, v ) || v.empty() || v[0] <= 0 ) writer->error("Cannot read sample interval '%s' from file '%s'.", dsDt.c_str(), f.name.c_str());
    f.dt_s = v[0];
    f.time_s = ( readNumbers( fid, dsTime, v ) && !v.empty() ) ? v[0] : -1.0;
    std::string dsChanNum = dsDt.substr( 0, dsDt.rfind('/') ) + "/channel_number";
    if( !readNumbers( fid, dsChanNum, f.chanNumber ) || f.chanNumber.size() != f.nChan ) f.chanNumber.clear();

    readNumbers( fid, grpMeta + "/REC_X", f.recX );
    readNumbers( fid, grpMeta + "/REC_Y", f.recY );
    readNumbers( fid, grpMeta + "/elevation", f.elev );
    readStrings( fid, grpMeta + "/sensor_id", f.sensorId );
    H5Fclose( fid );

    if( ifile > 0 ) {
      FileInfo const& f0 = vars->files[0];
      if( f.nSamp != f0.nSamp ) writer->error("File '%s' has %d samples, first file has %d.", f.name.c_str(), (int)f.nSamp, (int)f0.nSamp);
      if( fabs( f.dt_s - f0.dt_s ) > 1e-9 ) writer->error("File '%s' has a different sample interval.", f.name.c_str());
      if( f.nChan != f0.nChan ) writer->error("File '%s' has %d channels, first file has %d.", f.name.c_str(), (int)f.nChan, (int)f0.nChan);
    }
    writer->line( "  File %d: %s", (int)ifile+1, f.name.c_str() );
    writer->line( "    receivers: %d  channels: %d  samples: %d  dt: %g s  metadata: X/Y %s, elevation %s, sensor_id %s",
                  (int)f.nRec, (int)f.nChan, (int)f.nSamp, f.dt_s,
                  f.recX.size() == f.nRec ? "yes" : "no", f.elev.size() == f.nRec ? "yes" : "no",
                  f.sensorId.size() == f.nRec ? "yes" : "no" );
  }

  //--- Channel selection (user numbers are 1-based)
  int nChan = (int)vars->files[0].nChan;
  if( param->exists("channel") ) {
    int n = param->getNumValues( "channel" );
    for( int i = 0; i < n; i++ ) {
      int c;
      param->getInt( "channel", &c, i );
      if( c < 1 || c > nChan ) writer->error("Channel %d out of range (file has channels 1-%d).", c, nChan);
      vars->channels.push_back( c-1 );
    }
  }
  else {
    for( int c = 0; c < nChan; c++ ) vars->channels.push_back( c );
  }

  //--- Super header
  shdr->numSamples = (int)vars->files[0].nSamp;
  shdr->sampleInt  = (float)( vars->files[0].dt_s * 1000.0 );
  shdr->domain     = DOMAIN_XT;
  writer->line( "  Sample interval [ms]: %f", shdr->sampleInt );
  writer->line( "  Number of samples:    %d", shdr->numSamples );
  writer->line( "  Output order:         %s", vars->order == ORDER_CHANNEL ? "channel (all receivers of channel 1, then channel 2, ...)" : "receiver (all channels of receiver 1, then receiver 2, ...)" );

  //--- Trace headers
  vars->hdrId_trcno     = hdef->addStandardHeader( HDR_TRCNO.name );
  vars->hdrId_fileno    = hdef->addStandardHeader( HDR_FILENO.name );
  vars->hdrId_chan      = hdef->addStandardHeader( HDR_CHAN.name );
  vars->hdrId_rec_index = hdef->addStandardHeader( HDR_REC_INDEX.name );
  vars->hdrId_node      = hdef->addStandardHeader( HDR_NODE.name );
  vars->hdrId_rec_line  = hdef->addStandardHeader( HDR_REC_LINE.name );
  vars->hdrId_rcv       = hdef->addStandardHeader( HDR_RCV.name );
  vars->hdrId_rec_x     = hdef->addStandardHeader( HDR_REC_X.name );
  vars->hdrId_rec_y     = hdef->addStandardHeader( HDR_REC_Y.name );
  vars->hdrId_rec_z     = hdef->addStandardHeader( HDR_REC_Z.name );
  vars->hdrId_rec_elev  = hdef->addStandardHeader( HDR_REC_ELEV.name );
  vars->hdrId_rec_wdep  = hdef->addStandardHeader( HDR_REC_WDEP.name );
  vars->hdrId_time_samp1    = hdef->addStandardHeader( HDR_TIME_SAMP1.name );
  vars->hdrId_time_samp1_us = hdef->addStandardHeader( HDR_TIME_SAMP1_US.name );
  vars->hdrId_time_code = hdef->addStandardHeader( HDR_TIME_CODE.name );
  vars->hdrId_time_year = hdef->addStandardHeader( HDR_TIME_YEAR.name );
  vars->hdrId_time_day  = hdef->addStandardHeader( HDR_TIME_DAY.name );
  vars->hdrId_time_hour = hdef->addStandardHeader( HDR_TIME_HOUR.name );
  vars->hdrId_time_min  = hdef->addStandardHeader( HDR_TIME_MIN.name );
  vars->hdrId_time_sec  = hdef->addStandardHeader( HDR_TIME_SEC.name );
  vars->hdrId_sensor_id = hdef->addHeader( TYPE_STRING, "sensor_id", "Sensor ID string from HDF5 metadata", 64 );
}

//*************************************************************************************************
// Exec phase
//
void exec_mod_input_hdf5_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;

  if( vars->atEOF ) {
    traceGather->freeAllTraces();
    return;
  }

  //--- Open next file when needed
  if( vars->fileIndex < 0 || vars->dataId < 0 ) {
    vars->fileIndex += 1;
    if( vars->fileIndex >= (int)vars->files.size() ) {
      vars->atEOF = true;
      traceGather->freeAllTraces();
      return;
    }
    FileInfo const& f = vars->files[vars->fileIndex];
    vars->fileId = H5Fopen( f.name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
    if( vars->fileId < 0 ) writer->error("Cannot open HDF5 file '%s'.", f.name.c_str());
    vars->dataId  = H5Dopen2( vars->fileId, vars->dsData.c_str(), H5P_DEFAULT );
    vars->spaceId = H5Dget_space( vars->dataId );
    vars->recIndex = 0;
    vars->chanPos  = 0;
  }

  FileInfo const& f = vars->files[vars->fileIndex];
  hsize_t irec = vars->recIndex;
  int ichan    = vars->channels[vars->chanPos];

  //--- Read one trace: hyperslab [irec][ichan][0..nSamp-1]
  csTrace* trace = traceGather->trace(0);
  float* samples = trace->getTraceSamples();
  hsize_t start[3] = { irec, (hsize_t)ichan, 0 };
  hsize_t count[3] = { 1, 1, f.nSamp };
  H5Sselect_hyperslab( vars->spaceId, H5S_SELECT_SET, start, NULL, count, NULL );
  hid_t memspace = H5Screate_simple( 1, &count[2], NULL );
  herr_t status = H5Dread( vars->dataId, H5T_NATIVE_FLOAT, memspace, vars->spaceId, H5P_DEFAULT, samples );
  H5Sclose( memspace );
  if( status < 0 ) writer->error("Error reading receiver %d channel %d from file '%s'.", (int)irec+1, ichan+1, f.name.c_str());

  //--- Trace headers
  csTraceHeader* trcHdr = trace->getTraceHeader();
  vars->traceCounter += 1;
  trcHdr->setIntValue( vars->hdrId_trcno, vars->traceCounter );
  trcHdr->setIntValue( vars->hdrId_fileno, vars->fileIndex+1 );
  trcHdr->setIntValue( vars->hdrId_chan, f.chanNumber.empty() ? ichan+1 : (int)f.chanNumber[ichan] );
  trcHdr->setIntValue( vars->hdrId_rec_index, (int)irec+1 );
  if( f.recX.size() == f.nRec ) trcHdr->setDoubleValue( vars->hdrId_rec_x, f.recX[irec] );
  if( f.recY.size() == f.nRec ) trcHdr->setDoubleValue( vars->hdrId_rec_y, f.recY[irec] );
  if( f.elev.size() == f.nRec ) {
    trcHdr->setFloatValue( vars->hdrId_rec_z,    (float)f.elev[irec] );
    trcHdr->setFloatValue( vars->hdrId_rec_elev, (float)f.elev[irec] );
    trcHdr->setFloatValue( vars->hdrId_rec_wdep, (float)fabs( f.elev[irec] ) );
  }
  if( f.sensorId.size() == f.nRec ) {
    std::string const& id = f.sensorId[irec];
    trcHdr->setStringValue( vars->hdrId_sensor_id, id.substr( 0, 64 ) );
    std::string digits = ridDigits( id );
    if( !digits.empty() ) {
      trcHdr->setIntValue( vars->hdrId_rec_line, digitsToInt( digits, vars->ridLine1, vars->ridLine2 ) );
      trcHdr->setIntValue( vars->hdrId_rcv,      digitsToInt( digits, vars->ridSta1,  vars->ridSta2 ) );
      if( digits.size() <= 9 ) trcHdr->setIntValue( vars->hdrId_node, atoi( digits.c_str() ) );
    }
  }
  if( f.time_s >= 0 ) {
    double tsec = floor( f.time_s );
    time_t tt = (time_t)tsec;
    struct tm utc;
    gmtime_r( &tt, &utc );
    trcHdr->setIntValue( vars->hdrId_time_samp1,    (int)tsec );
    trcHdr->setIntValue( vars->hdrId_time_samp1_us, (int)floor( ( f.time_s - tsec ) * 1.0e6 + 0.5 ) );
    trcHdr->setIntValue( vars->hdrId_time_code, 4 );   // UTC
    trcHdr->setIntValue( vars->hdrId_time_year, utc.tm_year + 1900 );
    trcHdr->setIntValue( vars->hdrId_time_day,  utc.tm_yday + 1 );
    trcHdr->setIntValue( vars->hdrId_time_hour, utc.tm_hour );
    trcHdr->setIntValue( vars->hdrId_time_min,  utc.tm_min );
    trcHdr->setIntValue( vars->hdrId_time_sec,  utc.tm_sec );
  }

  //--- Advance to next trace
  bool fileDone = false;
  if( vars->order == ORDER_CHANNEL ) {
    vars->recIndex += 1;
    if( vars->recIndex >= f.nRec ) {
      vars->recIndex = 0;
      vars->chanPos += 1;
      if( vars->chanPos >= (int)vars->channels.size() ) fileDone = true;
    }
  }
  else {
    vars->chanPos += 1;
    if( vars->chanPos >= (int)vars->channels.size() ) {
      vars->chanPos = 0;
      vars->recIndex += 1;
      if( vars->recIndex >= f.nRec ) fileDone = true;
    }
  }
  if( fileDone ) {
    H5Sclose( vars->spaceId ); H5Dclose( vars->dataId ); H5Fclose( vars->fileId );
    vars->spaceId = vars->dataId = vars->fileId = -1;
  }
  (void)shdr;
}

//*************************************************************************************************
// Parameter definition (help: seaseis -m input_hdf5)
//
void params_mod_input_hdf5_( csParamDef* pdef ) {
  pdef->setModule( "INPUT_HDF5", "Input OBN data from HDF5 file(s)",
                   "Reads dataset [receiver][channel][sample]. Each receiver/channel pair becomes one trace. "
                   "Sets trace headers chan, rec_index, node, rec_line, rcv, rec_x, rec_y, rec_z, rec_elev, rec_wdep, sensor_id and time headers (UTC)" );

  pdef->addDoc("Default layout: /data [receiver][channel][sample] float, /header/dt [s], /header/time [s since 1970, UTC],");
  pdef->addDoc("/metadata/REC_X, /metadata/REC_Y, /metadata/elevation, /metadata/sensor_id (e.g. ..._RID116117011.P26774).");
  pdef->addDoc("From the digits after 'RID': rec_line = digits 1-4, rcv = digits 5-8, node = all digits (if <= 9 digits).");
  pdef->addDoc("Header chan = channel index (1, 2, ...), or the value in /header/channel_number if present (files written by OUTPUT_HDF5).");
  pdef->addDoc("Several files are read one after the other (all must have the same number of samples, channels and sample interval).");

  pdef->addParam( "filename", "Input HDF5 file name(s)", NUM_VALUES_VARIABLE, "Specify several file names on one line, or use several 'filename' lines" );
  pdef->addValue( "", VALTYPE_STRING, "Input file name" );

  pdef->addParam( "channel", "Channel(s) to read (1, 2, ...)", NUM_VALUES_VARIABLE, "Default: all channels" );
  pdef->addValue( "", VALTYPE_NUMBER, "Channel number (1 = first channel)" );

  pdef->addParam( "order", "Output trace order", NUM_VALUES_FIXED );
  pdef->addValue( "channel", VALTYPE_OPTION );
  pdef->addOption( "channel", "All receivers of the first selected channel, then all receivers of the next channel" );
  pdef->addOption( "receiver", "All selected channels of receiver 1, then receiver 2, ..." );

  pdef->addParam( "rid_line", "Digit positions of receiver line in the number after 'RID' in sensor_id", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "First digit (1 = first)" );
  pdef->addValue( "4", VALTYPE_NUMBER, "Last digit" );

  pdef->addParam( "rid_station", "Digit positions of receiver station in the number after 'RID' in sensor_id", NUM_VALUES_FIXED );
  pdef->addValue( "5", VALTYPE_NUMBER, "First digit (1 = first)" );
  pdef->addValue( "8", VALTYPE_NUMBER, "Last digit" );

  pdef->addParam( "dataset", "Name of the data set containing the samples", NUM_VALUES_FIXED );
  pdef->addValue( "/data", VALTYPE_STRING );

  pdef->addParam( "ds_dt", "Name of the data set containing the sample interval [s]", NUM_VALUES_FIXED );
  pdef->addValue( "/header/dt", VALTYPE_STRING );

  pdef->addParam( "ds_time", "Name of the data set containing the start time [s since 1970, UTC]", NUM_VALUES_FIXED );
  pdef->addValue( "/header/time", VALTYPE_STRING );

  pdef->addParam( "metadata", "Name of the group containing REC_X, REC_Y, elevation, sensor_id", NUM_VALUES_FIXED );
  pdef->addValue( "/metadata", VALTYPE_STRING );
}

//************************************************************************************************
// Start exec phase
//
bool start_exec_mod_input_hdf5_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

//************************************************************************************************
// Cleanup phase
//
void cleanup_mod_input_hdf5_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( vars->spaceId >= 0 ) H5Sclose( vars->spaceId );
  if( vars->dataId  >= 0 ) H5Dclose( vars->dataId );
  if( vars->fileId  >= 0 ) H5Fclose( vars->fileId );
  delete vars; vars = NULL;
}

extern "C" void _params_mod_input_hdf5_( csParamDef* pdef ) {
  params_mod_input_hdf5_( pdef );
}
extern "C" void _init_mod_input_hdf5_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_input_hdf5_( param, env, writer );
}
extern "C" bool _start_exec_mod_input_hdf5_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_input_hdf5_( env, writer );
}
extern "C" void _exec_mod_input_hdf5_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_input_hdf5_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_input_hdf5_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_input_hdf5_( env, writer );
}
