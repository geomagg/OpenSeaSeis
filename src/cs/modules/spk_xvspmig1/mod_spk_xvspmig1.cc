/* Copyright (c) 2026 M. A. Gallotti Guimaraes (COPPE-UFRJ).      */
/* Distributed under the OpenSeaSeis license (see LICENSE).       */
/* Module SPK_XVSPMIG1: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_XVSPMIG1
 *
 * 2D prestack depth migration, record by record, from SEISPAK program XVSPMIG1
 * (sgivspmig1 / SHOTMIG) for the UP-GOING field of OBN receiver gathers: by
 * reciprocity the traces are sources at the surface and the node is a point source
 * at depth GZ. The recorded field is continued down from the surface; the source
 * field starts two steps below the node; both with multi-velocity PSPI. Imaging by
 * correlation below the node (no mirror), amplitude restoration, angle mute starting
 * at the node depth (20-sample ramp) and stack of all records.
 */

extern "C" {
  void spkvsp1_( float* d, int* nt, int* nntr, int* nrecs, float* sr, float* afit, float* alph, float* gz, float* spjn,
                int* nnz, int* itzr, float* ddz, float* dx, float* apx, float* frq1, float* frq2, float* oty,
                float* amute, int* ltm, int* lev1, float* stk, int* nstk, float* outr, int* ioutr, int* info, int* ierr,
                int* ifoc, float* rfoc, float* tdat, float* otau, float* ovc, float* ovavg, float* picks );
  void spkvxt2_( float* vel, float* dx, int* nxv, int* nyv, int* np, float* afit, float* ee, float* x, float* bfit, float* xx );
}

namespace mod_spk_xvspmig1 {
  struct VariableStruct {
    int   ntIn, nnz, itzr, ltm, lev1;
    float sr, ddz, dx, apx, frq1, frq2, oty, amute;
    bool  geomHeaders;
    std::string recHeader, hdrTraceX, hdrTraceY, hdrSrcX, hdrSrcY, hdrSrcZ;
    float alph1, gam1, alph2, alph3, zzy;
    bool  outputRecords;
    std::vector<float> afit;
    int   nxv, np;
    csTraceGather* gather;
    // coherence (velocity focussing) analysis of xvspmig1: FOR LOR RNI WINDOW VMIN VMAX EPS
    bool  focus; int ifor, ilor, irni; float window, vmin, vmax, seps;
    std::string focusOut;
  };
  // SEASEIS: small SEG-Y writer (rev 1, IEEE big-endian) for the coherence analysis outputs
  static void putBI( unsigned char* p, int v ) { p[0]=(v>>24)&255; p[1]=(v>>16)&255; p[2]=(v>>8)&255; p[3]=v&255; }
  static void putBS( unsigned char* p, int v ) { p[0]=(v>>8)&255; p[1]=v&255; }
  struct SegyOut {
    FILE* fp; int ns; int count;
    SegyOut() : fp(NULL), ns(0), count(0) {}
    bool open( std::string const& fn, int nsamp, int sampleUs ) {
      fp = fopen( fn.c_str(), "wb" ); ns = nsamp; count = 0;
      if( fp == NULL ) return false;
      unsigned char ebc[3200]; memset( ebc, 0x40, 3200 ); fwrite( ebc, 1, 3200, fp );
      unsigned char bh[400]; memset( bh, 0, 400 );
      putBI( bh, 1 ); putBS( bh+16, sampleUs ); putBS( bh+20, nsamp ); putBS( bh+24, 5 ); putBS( bh+54, 1 );
      putBS( bh+300, 0x0100 ); putBS( bh+302, 1 );
      fwrite( bh, 1, 400, fp );
      return true;
    }
    void trace( float const* d, int sampleUs, int fldr, int trc, int key, int cdp ) {
      unsigned char th[240]; memset( th, 0, 240 );
      count++;
      putBI( th, count ); putBI( th+4, count ); putBI( th+8, fldr ); putBI( th+12, trc ); putBI( th+20, cdp );
      putBS( th+28, 1 ); putBI( th+36, key ); putBS( th+70, 1 ); putBS( th+114, ns ); putBS( th+116, sampleUs );
      putBI( th+188, 1 ); putBI( th+192, cdp );
      fwrite( th, 1, 240, fp );
      std::vector<unsigned char> b( 4*ns );
      for( int i = 0; i < ns; i++ ) { unsigned int u; memcpy( &u, &d[i], 4 ); putBI( &b[4*i], (int)u ); }
      fwrite( &b[0], 1, b.size(), fp );
    }
    void close() { if( fp != NULL ) fclose( fp ); fp = NULL; }
  };

  static int const BASE = 12;
  static inline float getF( char const* p, bool sw ) {
    char b[4]; memcpy(b,p,4); if(sw){ char c=b[0];b[0]=b[3];b[3]=c; c=b[1];b[1]=b[2];b[2]=c; }
    float v; memcpy(&v,b,4); return v;
  }
  static inline int getI( char const* p, bool sw ) {
    char b[4]; memcpy(b,p,4); if(sw){ char c=b[0];b[0]=b[3];b[3]=c; c=b[1];b[1]=b[2];b[2]=c; }
    int v; memcpy(&v,b,4); return v;
  }
  static inline short getS( char const* p, bool sw ) {
    char b[2]; memcpy(b,p,2); if(sw){ char c=b[0];b[0]=b[1];b[1]=c; }
    short v; memcpy(&v,b,2); return v;
  }
  std::string readSeispak( std::string const& fname, int* hmtOut, int* hmrOut, int* ntOut, float* smpOut,
                           std::vector<float>* data ) {
    FILE* fp = fopen( fname.c_str(), "rb" );
    if( fp == NULL ) return "cannot open file";
    std::vector<char> raw;
    char buf[65536];
    size_t n;
    while( (n = fread(buf,1,sizeof(buf),fp)) > 0 ) raw.insert( raw.end(), buf, buf+n );
    fclose(fp);
    if( raw.size() < 200 || strncmp(&raw[0],"ICE",3) != 0 ) return "not a SEISPAK file (no 'ICE' signature)";
    char const* d = &raw[BASE];
    bool sw = false;
    int lh = getI(d+4,false);
    if( !(lh >= 0 && lh < 10000) ) { sw = true; lh = getI(d+4,true); }
    int ifw0 = getI(d,sw), nt = getI(d+8,sw), nbyte = getI(d+32,sw);
    float smp = getF(d+12,sw), hmt = getF(d+16,sw), hmr = getF(d+20,sw);
    if( nbyte == 0 ) nbyte = 4;
    int ihmt = (int)(hmt+0.5f), ihmr = (int)(hmr+0.5f);
    int ntr = ihmt*ihmr;
    if( nt <= 0 || ntr <= 0 || ifw0 < 21 ) return "implausible header directory";
    long long rec;
    if( nbyte == 4 ) rec = 4LL*(lh+nt);
    else if( nbyte == 2 ) rec = 4LL*(lh + nt/2 + ((nt%2==0)?1:2));
    else return "only NBYTE 2 and 4 supported for velocity files";
    long long first = BASE + 4LL*(ifw0-1);
    if( first + rec*ntr > (long long)raw.size() ) return "file truncated";
    data->assign( (size_t)ntr*nt, 0.0f );
    for( int k = 0; k < ntr; k++ ) {
      char const* p = &raw[first + rec*k] + 4*lh;
      float* out = &(*data)[(size_t)k*nt];
      if( nbyte == 4 ) {
        for( int i = 0; i < nt; i++ ) out[i] = getF(p+4*i,sw);
      }
      else {
        float fac = getF(p,sw);
        char const* s = p + 4 + ((nt%2) ? 2 : 0);
        for( int i = 0; i < nt; i++ ) out[i] = (fac > 1e-30f) ? getS(s+2*i,sw)/fac : 0.0f;
      }
    }
    *hmtOut = ihmt; *hmrOut = ihmr; *ntOut = nt; *smpOut = smp;
    return "";
  }
  // SEISPAK onl2dvsvels: v(z) of each DPN (sample interval dzv) -> NZ interval velocities at depth
  // step DZ, preserving the vertical travel time (t(z) integrated with step dzv, linear
  // interpolation of t at j*DZ (int3), v = DZ/dt). Pairs (v_j, j*DZ), as written to DVSVELS.
  void onl2dvsvels( std::vector<float> const& vin, int ntr, int nzv, float dzv, int nz, float dz, std::vector<float>* vout ) {
    int nnz = (int)lrint( (nz+1)*(double)dz/dzv ) + 2;
    std::vector<double> y( nnz+2 ), f( nz+2 );
    double r = (double)dz/dzv;
    int m = nz+1;
    int m1 = std::min( m, 1 + (int)( (nnz-1)/r ) );
    vout->assign( (size_t)ntr*nz, 0.0f );
    for( int ix = 0; ix < ntr; ix++ ) {
      float const* v = &vin[(size_t)ix*nzv];
      double t = 0.0;
      for( int j = 1; j <= nnz+1; j++ ) {
        int iv = std::min( nzv, j );
        y[j] = t;
        t += dzv/v[iv-1];
      }
      for( int j = 1; j <= m1; j++ ) {
        double x = 1.0 + r*(j-1);
        int ii = (int)x;
        double b = x - ii, a = 1.0 - b;
        f[j] = a*y[ii] + b*y[ii+1];
      }
      for( int j = m1+1; j <= m; j++ ) f[j] = 0.0;
      double t1 = 0.0;
      for( int j = 1; j <= nz; j++ ) {
        double t2 = f[j+1];
        (*vout)[(size_t)ix*nz + j-1] = (float)( dz/(t2-t1) );
        t1 = t2;
      }
    }
  }
  double hval( csTrace* trc, int id ) {
    return ( id < 0 ) ? 0.0 : trc->getTraceHeader()->doubleValue( id );
  }
}
using namespace mod_spk_xvspmig1;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_spk_xvspmig1_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  std::string text;
  vars->apx = 12.0f; vars->itzr = 5; vars->frq1 = 5.0f; vars->frq2 = 60.0f; vars->oty = 0.0f;
  vars->amute = 30.0f; vars->ltm = 0; vars->lev1 = 0; vars->nnz = 0; vars->ddz = 0.0f; vars->dx = 0.0f;
  vars->geomHeaders = true;
  vars->recHeader = "";
  vars->hdrTraceX = "sou_x"; vars->hdrTraceY = "sou_y";
  vars->hdrSrcX = "rec_x";   vars->hdrSrcY = "rec_y"; vars->hdrSrcZ = "rec_z";
  vars->alph1 = 0.0f; vars->gam1 = 0.0f; vars->alph2 = 0.0f; vars->alph3 = 0.0f; vars->zzy = 0.0f;
  vars->outputRecords = false;
  vars->nxv = vars->np = 0;
  vars->gather = NULL;

  if( param->exists("apx") )   param->getFloat( "apx", &vars->apx );
  if( param->exists("itzr") )  param->getInt( "itzr", &vars->itzr );
  if( param->exists("frq1") )  param->getFloat( "frq1", &vars->frq1 );
  if( param->exists("frq2") )  param->getFloat( "frq2", &vars->frq2 );
  if( param->exists("dx") )    param->getFloat( "dx", &vars->dx );
  if( param->exists("ddz") )   param->getFloat( "ddz", &vars->ddz );
  if( param->exists("nnz") )   param->getInt( "nnz", &vars->nnz );
  if( param->exists("oty") )   param->getFloat( "oty", &vars->oty );
  if( param->exists("mute_angle") ) param->getFloat( "mute_angle", &vars->amute );
  if( param->exists("ltm") )   param->getInt( "ltm", &vars->ltm );
  if( param->exists("lev") )   param->getInt( "lev", &vars->lev1 );
  if( param->exists("record_header") ) param->getString( "record_header", &vars->recHeader );
  if( param->exists("geometry") ) {
    param->getString( "geometry", &text );
    if( !text.compare("parameters") ) vars->geomHeaders = false;
    else if( text.compare("headers") ) writer->error("Unknown option for 'geometry': %s", text.c_str());
  }
  if( param->exists("hdr_trace_x") )  param->getString( "hdr_trace_x", &vars->hdrTraceX );
  if( param->exists("hdr_trace_y") )  param->getString( "hdr_trace_y", &vars->hdrTraceY );
  if( param->exists("hdr_source_x") ) param->getString( "hdr_source_x", &vars->hdrSrcX );
  if( param->exists("hdr_source_y") ) param->getString( "hdr_source_y", &vars->hdrSrcY );
  if( param->exists("hdr_source_z") ) param->getString( "hdr_source_z", &vars->hdrSrcZ );
  if( param->exists("alph1") ) param->getFloat( "alph1", &vars->alph1 );
  if( param->exists("gam1") )  param->getFloat( "gam1", &vars->gam1 );
  if( param->exists("alph2") ) param->getFloat( "alph2", &vars->alph2 );
  if( param->exists("alph3") ) param->getFloat( "alph3", &vars->alph3 );
  if( param->exists("zzy") )   param->getFloat( "zzy", &vars->zzy );
  vars->focus = false; vars->ifor = 0; vars->ilor = -1; vars->irni = 20; vars->window = 1000.0f;
  vars->vmin = 5000.0f; vars->vmax = 10000.0f; vars->seps = 0.1f; vars->focusOut = "FOCUS";
  if( param->exists("focus_x") ) {
    vars->focus = true;
    param->getInt( "focus_x", &vars->ifor, 0 ); param->getInt( "focus_x", &vars->ilor, 1 );
    if( param->getNumValues("focus_x") > 2 ) param->getInt( "focus_x", &vars->irni, 2 );
    if( vars->irni < 1 || vars->ilor < vars->ifor ) writer->error("focus_x must be FOR LOR RNI with LOR >= FOR and RNI >= 1");
    if( param->exists("focus_window") ) param->getFloat( "focus_window", &vars->window );
    if( param->exists("focus_vel") ) { param->getFloat( "focus_vel", &vars->vmin, 0 ); param->getFloat( "focus_vel", &vars->vmax, 1 ); }
    if( vars->vmax <= vars->vmin ) writer->error("focus_vel must be VMIN VMAX with VMIN < VMAX");
    if( param->exists("focus_eps") ) param->getFloat( "focus_eps", &vars->seps );
    if( param->exists("focus_out") ) param->getString( "focus_out", &vars->focusOut );
    writer->line("  Coherence analysis (xvspmig1): DPN %d-%d step %d, window %g ms, velocities %g-%g, EPS %g; output '%s_tau.sgy' (OFFLINE4), '%s_vcoh.sgy' (OFFLINE5), '%s_vavg.sgy' (VAVG), '%s_vtps.txt'",
                 vars->ifor, vars->ilor, vars->irni, vars->window, vars->vmin, vars->vmax, vars->seps,
                 vars->focusOut.c_str(), vars->focusOut.c_str(), vars->focusOut.c_str(), vars->focusOut.c_str());
  }
  if( param->exists("output") ) {
    param->getString( "output", &text );
    if( !text.compare("records") ) vars->outputRecords = true;
    else if( text.compare("stack") ) writer->error("Unknown option for 'output': %s", text.c_str());
  }

  if( vars->apx < 11.0f ) writer->error("APX = %g not supported: XVSPMIG1 is ported for multi-velocity PSPI only, APX >= 11 (nvels = APX-10)", vars->apx);
  if( vars->dx <= 0.0f )  writer->error("Parameter 'dx' (trace spacing within a record) is required");
  if( vars->ddz <= 0.0f || vars->nnz <= 0 ) writer->error("Parameters 'ddz' and 'nnz' (depth sampling of the image) are required");
  if( vars->itzr < 1 ) writer->error("ITZR must be >= 1");
  if( vars->frq2 <= vars->frq1 ) writer->error("FRQ2 must be larger than FRQ1");
  if( vars->ltm < 0 ) writer->error("LTM must be >= 0");
  if( !vars->geomHeaders && vars->alph1 <= 0.0f ) writer->error("geometry parameters: ALPH1 (trace position of the source, 1 = first trace) is required");

  if( vars->recHeader.empty() ) {
    if( hdef->headerExists("ffid") ) vars->recHeader = "ffid";
    else if( hdef->headerExists("spk_rec") ) vars->recHeader = "spk_rec";
    else writer->error("Cannot find records: no header 'ffid' or 'spk_rec'. Use 'record_header'.");
  }
  else if( !hdef->headerExists(vars->recHeader) ) writer->error("Record header '%s' does not exist", vars->recHeader.c_str());
  if( vars->geomHeaders ) {
    if( !hdef->headerExists(vars->hdrTraceX) ) writer->error("Trace position header '%s' does not exist (use hdr_trace_x or geometry parameters)", vars->hdrTraceX.c_str());
    if( !hdef->headerExists(vars->hdrSrcX) )   writer->error("Source position header '%s' does not exist (use hdr_source_x or geometry parameters)", vars->hdrSrcX.c_str());
    if( !hdef->headerExists(vars->hdrSrcZ) )   writer->error("Source depth header '%s' does not exist (use hdr_source_z or geometry parameters)", vars->hdrSrcZ.c_str());
  }

  vars->ntIn = shdr->numSamples;
  vars->sr   = shdr->sampleInt;
  float fnyq = 500.0f / vars->sr;
  if( vars->frq2 > fnyq ) writer->warning("FRQ2 (%f) above Nyquist (%f Hz): limited to Nyquist", vars->frq2, fnyq);

  //---------------------------------------------
  // Velocity -> dense functions: afit(1..3) = nxv, 1, np; per location iy, ix, pairs (v, z)
  bool hasVxt = param->exists("vxt_file"), hasVel = param->exists("vel_file");
  if( hasVxt == hasVel ) writer->error("Specify one velocity source: 'vxt_file' (SEISPAK VXT list) or 'vel_file' (SEISPAK ONL v(z) file)");
  if( hasVxt ) {
    std::string fname;
    param->getString( "vxt_file", &fname );
    float yLine = -1.0f;
    if( param->exists("vxt_line") ) param->getFloat( "vxt_line", &yLine );
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
      if( ig < 0 ) writer->error("VXT file: line y=%g not found", yLine);
    }
    std::vector<float> list;
    list.push_back( v[0] );
    list.insert( list.end(), lead.begin(), lead.end() );
    list.push_back( 1.0f );
    list.insert( list.end(), bodies[ig].begin(), bodies[ig].end() );
    list.push_back( 0.0f );
    size_t m = list.size() + 60000;
    list.resize( m, 0.0f );
    std::vector<float> ee( m ), x( 20*m + 100000 ), xx( 10001 );
    size_t nafit = 3 + (size_t)10001*203;
    vars->afit.assign( nafit, 0.0f );
    std::vector<float> bfit( nafit );
    int nxv = 0, nyv = 0, np = 0;
    float ddd = vars->dx;
    spkvxt2_( &list[0], &ddd, &nxv, &nyv, &np, &vars->afit[0], &ee[0], &x[0], &bfit[0], &xx[0] );
    if( nxv <= 0 || np <= 0 || nxv > 10000 ) writer->error("VXT file '%s': could not build velocity functions (nx=%d, events=%d)", fname.c_str(), nxv, np);
    vars->nxv = nxv; vars->np = np;
    writer->line("  VXT file: %s  (%d y-lines in file, using y = %g), %d horizons, %d DPN", fname.c_str(), (int)ys.size(), ys[ig], np, nxv);
  }
  else {
    std::string fname;
    param->getString( "vel_file", &fname );
    int hmt, hmr, ntv; float dzv;
    std::vector<float> vdat;
    std::string err = readSeispak( fname, &hmt, &hmr, &ntv, &dzv, &vdat );
    if( !err.empty() ) {
      std::string msg = "Velocity file '" + fname + "': " + err;
      writer->error( "%s", msg.c_str() );
    }
    if( param->exists("vel_dz") ) {                // onl2dvsvels NX NY NZ DZ
      float dz = 0.0f; int nz = 0;
      param->getFloat( "vel_dz", &dz );
      if( dz <= 0.0f ) writer->error("vel_dz must be > 0");
      if( dzv <= 0.0f ) writer->error("Velocity file '%s': sample interval %g in its header is not a depth step", fname.c_str(), dzv);
      if( param->exists("vel_nz") ) param->getInt( "vel_nz", &nz );
      else nz = (int)( ntv*(double)dzv/dz ) - 1;
      if( nz < 1 ) writer->error("vel_nz must be >= 1");
      std::vector<float> vv;
      onl2dvsvels( vdat, hmt*hmr, ntv, dzv, nz, dz, &vv );
      writer->line("  onl2dvsvels: %d samples at %g -> %d samples at %g (vertical time preserved)", ntv, dzv, nz, dz);
      vdat.swap( vv ); ntv = nz; dzv = dz;
    }
    int nxv = hmt*hmr, np = ntv, nw = 2*np+3;
    vars->afit.assign( 3 + (size_t)nxv*nw + 1000, 0.0f );
    for( int ix = 1; ix <= nxv; ix++ ) {
      float* xx = &vars->afit[3 + (size_t)(ix-1)*nw];
      float const* tr = &vdat[(size_t)(ix-1)*ntv];
      xx[0] = 1.0f; xx[1] = (float)ix;
      for( int j = 0; j < np; j++ ) { xx[2+2*j] = tr[j]; xx[3+2*j] = (j+1)*dzv; }
    }
    vars->afit[0] = (float)nxv; vars->afit[1] = 1.0f; vars->afit[2] = (float)np;
    vars->nxv = nxv; vars->np = np;
    writer->line("  Velocity file (ONL): %s  (%d DPN, %d samples, depth step %f)", fname.c_str(), nxv, np, dzv);
  }

  //---------------------------------------------
  shdr->numSamples = vars->nnz + vars->ltm;
  shdr->sampleInt  = vars->ddz;
  shdr->domain     = DOMAIN_XD;
  vars->gather = new csTraceGather( hdef );
  if( !vars->outputRecords ) {
    if( !hdef->headerExists(HDR_TRCNO.name) ) hdef->addStandardHeader( HDR_TRCNO.name );
    if( !hdef->headerExists(HDR_CMP.name) )   hdef->addStandardHeader( HDR_CMP.name );
    if( !hdef->headerExists(HDR_CMP_X.name) ) hdef->addStandardHeader( HDR_CMP_X.name );
  }

  writer->line("  APX %g (multi-velocity PSPI, %d velocities)  ITZR %d  FRQ1/FRQ2 %g %g Hz", vars->apx, (int)lrintf(vars->apx-10.0f), vars->itzr, vars->frq1, vars->frq2);
  writer->line("  DX %g  DDZ %g  NNZ %d  LTM %d  OTY %g  mute angle %g deg", vars->dx, vars->ddz, vars->nnz, vars->ltm, vars->oty, vars->amute);
  writer->line("  Records from header '%s'", vars->recHeader.c_str());
  if( vars->geomHeaders ) writer->line("  Geometry from headers: traces '%s', source (node) '%s' at depth '%s'", vars->hdrTraceX.c_str(), vars->hdrSrcX.c_str(), vars->hdrSrcZ.c_str());
  else writer->line("  Geometry from parameters: ALPH1 %g (+GAM1 %g per record), depth ALPH2 %g (+ALPH3 %g), ZZY %g", vars->alph1, vars->gam1, vars->alph2, vars->alph3, vars->zzy);
  writer->line("  Output: %s", vars->outputRecords ? "migrated records (input headers)" : "stacked depth section");
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_spk_xvspmig1_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csTraceHeaderDef const* hdef = env->headerDef;

  if( !edef->isLastCall() ) {
    traceGather->moveTracesTo( 0, traceGather->numTraces(), vars->gather );
    edef->setTracesAreWaiting();
    return;
  }
  int ntr = vars->gather->numTraces();
  if( ntr == 0 ) return;

  //--------------------------------------------------------------
  // Records
  int idRec = hdef->headerIndex( vars->recHeader );
  std::vector<int> recStart;
  double prev = 0.0;
  for( int k = 0; k < ntr; k++ ) {
    double v = hval( vars->gather->trace(k), idRec );
    if( k == 0 || v != prev ) recStart.push_back( k );
    prev = v;
  }
  int nrecs = (int)recStart.size();
  recStart.push_back( ntr );
  int nntr = recStart[1] - recStart[0];
  for( int i = 1; i < nrecs; i++ ) {
    if( recStart[i+1] - recStart[i] != nntr ) {
      char msg[256];
      snprintf( msg, sizeof(msg), "Record %d has %d traces, first record has %d. All records must have the same number of traces.", i+1, recStart[i+1]-recStart[i], nntr );
      writer->error( "%s", msg );
    }
  }
  if( nntr < 4 ) writer->error("Records must have at least 4 traces");
  writer->line("  Data: %d records x %d traces, %d samples", nrecs, nntr, vars->ntIn);

  //--------------------------------------------------------------
  // Geometry per record: ALPH1 (source trace position), GZ (source depth), SPJN (shift of trace 1)
  std::vector<float> alph( nrecs ), gz( nrecs ), spjn( nrecs );
  double x0 = 0.0;
  if( vars->geomHeaders ) {
    int idTx = hdef->headerIndex( vars->hdrTraceX ), idTy = hdef->headerExists( vars->hdrTraceY ) ? hdef->headerIndex( vars->hdrTraceY ) : -1;
    int idSx = hdef->headerIndex( vars->hdrSrcX ),   idSy = hdef->headerExists( vars->hdrSrcY ) ? hdef->headerIndex( vars->hdrSrcY ) : -1;
    int idSz = hdef->headerIndex( vars->hdrSrcZ );
    // line direction from first to last trace of record 1
    csTrace* a = vars->gather->trace(0);
    csTrace* b = vars->gather->trace(nntr-1);
    double ux = hval(b,idTx) - hval(a,idTx), uy = hval(b,idTy) - hval(a,idTy);
    double len = sqrt( ux*ux + uy*uy );
    if( len <= 0.0 ) writer->error("Header '%s' is the same on all traces of record 1: cannot find the line direction. Use geometry parameters.", vars->hdrTraceX.c_str());
    ux /= len; uy /= len;
    double spacing = len / (nntr-1);
    if( fabs( spacing - vars->dx ) > 0.05*vars->dx ) writer->warning("Mean trace spacing in record 1 from headers is %g, but DX = %g", spacing, vars->dx);
    std::vector<double> p1( nrecs );
    for( int i = 0; i < nrecs; i++ ) {
      csTrace* t1 = vars->gather->trace( recStart[i] );
      csTrace* tn = vars->gather->trace( recStart[i+1]-1 );
      double dirx = hval(tn,idTx) - hval(t1,idTx), diry = hval(tn,idTy) - hval(t1,idTy);
      if( dirx*ux + diry*uy <= 0.0 ) writer->error("Traces of record %d are not ordered in the same direction as record 1 (sort each record by trace position)", i+1);
      p1[i] = hval(t1,idTx)*ux + hval(t1,idTy)*uy;
      double ps = hval(t1,idSx)*ux + hval(t1,idSy)*uy;
      alph[i] = (float)( (ps - p1[i]) / vars->dx + 1.0 );
      gz[i]   = (float)fabs( hval(t1,idSz) );
    }
    double pmin = *std::min_element( p1.begin(), p1.end() );
    for( int i = 0; i < nrecs; i++ ) spjn[i] = (float)( (p1[i] - pmin) / vars->dx );
    x0 = pmin;
  }
  else {
    for( int i = 0; i < nrecs; i++ ) {
      alph[i] = vars->alph1 + i*vars->gam1;
      gz[i]   = vars->alph2 + i*vars->alph3;
      spjn[i] = i*vars->zzy;
    }
  }
  for( int i = 0; i < nrecs; i++ ) {
    if( alph[i] < -15.0f || alph[i] > nntr + 16.0f ) writer->warning("Record %d: source at trace %.1f, outside the traces (1..%d): its image will be empty", i+1, alph[i], nntr);
  }
  writer->line("  Record 1: source at trace %.2f, depth %g;  record %d: source at trace %.2f, depth %g, trace 1 shifted %.2f traces",
               alph[0], gz[0], nrecs, alph[nrecs-1], gz[nrecs-1], spjn[nrecs-1]);

  //--------------------------------------------------------------
  int nt = vars->ntIn;
  std::vector<float> d( (size_t)ntr*nt );
  for( int k = 0; k < ntr; k++ ) memcpy( &d[(size_t)k*nt], vars->gather->trace(k)->getTraceSamples(), nt*sizeof(float) );
  float spmax = *std::max_element( spjn.begin(), spjn.end() );
  int nnzt = vars->nnz + vars->ltm;
  int nstk = nntr + (int)(spmax + 0.5f);
  std::vector<float> stk( (size_t)nnzt*nstk, 0.0f );
  int ioutr = vars->outputRecords ? 1 : 0;
  std::vector<float> outr( ioutr ? (size_t)nnzt*ntr : 1, 0.0f );
  int info[8] = {0,0,0,0,0,0,0,0};
  int ierr = 0;
  // coherence analysis arrays (TDAT(jnt,nz2,ngsn) and the outputs)
  int jnt = vars->focus ? (int)(vars->window/vars->sr) + 1 : 1;
  int ngsn = vars->focus ? (vars->ilor - vars->ifor)/vars->irni + 1 : 1;
  int nz2 = (vars->nnz - 1)/vars->itzr + 1 + vars->ltm/vars->itzr;
  int iisr = std::max( 1, (int)(8.0f/vars->sr + 0.1f) ), jnt1 = (jnt-1)/iisr + 1;
  int ifoc[4] = { vars->focus ? vars->ifor : 0, vars->focus ? vars->ilor : -1, vars->irni, jnt };
  float rfoc[3] = { vars->vmin, vars->vmax, vars->seps };
  std::vector<float> tdat( vars->focus ? (size_t)jnt*nz2*ngsn + 10 : 1 ), otau( vars->focus ? (size_t)nt*jnt1*ngsn : 1 ),
                     ovc( otau.size() ), ovavg( vars->focus ? (size_t)nt*ngsn : 1 ), picks( vars->focus ? (size_t)201*ngsn : 1, 0.0f );
  float frq2 = vars->frq2;
  writer->line("  Migrating...");
  spkvsp1_( &d[0], &nt, &nntr, &nrecs, &vars->sr, &vars->afit[0], &alph[0], &gz[0], &spjn[0], &vars->nnz, &vars->itzr,
           &vars->ddz, &vars->dx, &vars->apx, &vars->frq1, &frq2, &vars->oty, &vars->amute, &vars->ltm, &vars->lev1,
           &stk[0], &nstk, &outr[0], &ioutr, info, &ierr,
           ifoc, rfoc, &tdat[0], &otau[0], &ovc[0], &ovavg[0], &picks[0] );
  writer->line("  nx (FFT) %d  nt1 %d  frequency samples %d-%d  depth steps %d  velocity grid %d  memory %d MB",
               info[0], info[1], info[2], info[3], info[4], info[5], info[6]);
  if( ierr == 1 ) writer->error("No frequencies to migrate between FRQ1 and FRQ2");
  if( ierr == 2 ) writer->error("Not enough memory for migration (%d MB)", info[6]);
  if( vars->focus ) {
    std::string b = vars->focusOut;
    int sus = (int)lrint( vars->sr*1000.0 ), mjnt = jnt/2 + 1;
    float dv = (vars->vmax - vars->vmin)/(jnt - 1);
    SegyOut ft, fv, fa;
    if( !ft.open( b + "_tau.sgy", nt, sus ) || !fv.open( b + "_vcoh.sgy", nt, sus ) || !fa.open( b + "_vavg.sgy", nt, sus ) )
      writer->error("Cannot write the coherence analysis files '%s_*.sgy'", b.c_str());
    FILE* fvt = fopen( (b + "_vtps.txt").c_str(), "w" );
    if( fvt != NULL ) { fprintf( fvt, " &rmspars\n vtp = \n" ); }
    for( int i = 1; i <= ngsn; i++ ) {
      int dpn = vars->ifor + vars->irni*(i-1);
      for( int k = 1; k <= jnt1; k++ ) {
        int j = (k-1)*iisr + 1;
        size_t o = ((size_t)(i-1)*jnt1 + (k-1))*nt;
        ft.trace( &otau[o], sus, i, k, (int)lrint((j - mjnt)*vars->sr), dpn );
        fv.trace( &ovc[o], sus, i, k, (int)lrint(vars->vmin + (j-1)*dv), dpn );
      }
      fa.trace( &ovavg[(size_t)(i-1)*nt], sus, i, 1, 0, dpn );
      float const* pk = &picks[(size_t)(i-1)*201];
      int npk = (int)lrint( pk[0] );
      if( fvt != NULL && npk > 0 ) {
        fprintf( fvt, "%d\n", dpn );
        for( int j = 1; j <= npk; j++ ) fprintf( fvt, "%d %d\n", (int)pk[2*j-1], (int)pk[2*j] );
        fprintf( fvt, "0\n" );
      }
    }
    if( fvt != NULL ) { fprintf( fvt, " &end\n" ); fclose( fvt ); }
    ft.close(); fv.close(); fa.close();
    writer->line("  Coherence analysis: %d points x %d traces written to '%s_tau.sgy', '%s_vcoh.sgy'; VAVG '%s_vavg.sgy'; picks '%s_vtps.txt'",
                 ngsn, jnt1, b.c_str(), b.c_str(), b.c_str(), b.c_str());
  }
  if( ierr == 3 ) writer->error("Node depth (%d depth steps) is too close to the image bottom: increase NNZ", info[7]);

  //--------------------------------------------------------------
  if( vars->outputRecords ) {
    traceGather->createTraces( 0, ntr, hdef, nnzt );
    for( int k = 0; k < ntr; k++ ) {
      csTrace* trc = traceGather->trace(k);
      trc->getTraceHeader()->copyFrom( vars->gather->trace(k)->getTraceHeader() );
      memcpy( trc->getTraceSamples(), &outr[(size_t)k*nnzt], nnzt*sizeof(float) );
    }
  }
  else {
    // last trace with data
    int nout = nstk;
    while( nout > 1 ) {
      bool any = false;
      for( int j = 0; j < nnzt && !any; j++ ) if( stk[(size_t)(nout-1)*nnzt + j] != 0.0f ) any = true;
      if( any ) break;
      nout--;
    }
    traceGather->createTraces( 0, nout, hdef, nnzt );
    int idTrc = hdef->headerIndex( HDR_TRCNO.name );
    int idCmp = hdef->headerIndex( HDR_CMP.name );
    int idCmpX = hdef->headerIndex( HDR_CMP_X.name );
    for( int k = 0; k < nout; k++ ) {
      csTrace* trc = traceGather->trace(k);
      csTraceHeader* th = trc->getTraceHeader();
      th->setIntValue( idTrc, k+1 );
      th->setIntValue( idCmp, k+1 );
      th->setDoubleValue( idCmpX, x0 + k*vars->dx );
      memcpy( trc->getTraceSamples(), &stk[(size_t)k*nnzt], nnzt*sizeof(float) );
    }
    writer->line("  Output: %d stacked traces", nout);
    if( vars->nxv < nout ) writer->warning("Velocity model has %d DPN, the section has %d traces: the last DPN is repeated", vars->nxv, nout);
  }
  vars->gather->freeAllTraces();
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_spk_xvspmig1_( csParamDef* pdef ) {
  pdef->setModule( "SPK_XVSPMIG1", "2D prestack depth migration of the up-going field of OBN receiver gathers (SEISPAK XVSPMIG1 / SHOTMIG)",
                   "Each record (receiver gather of one node; traces = shots at the surface) is migrated by reciprocity with a point "
                   "source at the node (trace position ALPH1, depth GZ). The recorded up-going field (e.g. PZ separated) is continued down "
                   "from the surface, the source field starts below the node, both with multi-velocity PSPI; imaging by correlation below "
                   "the node (no mirror image), amplitude restoration, angle mute from the node depth down, and stack of all records. "
                   "For the down-going field (mirror imaging) use SPK_XVSPMIG2. "
                   "Input: records with the same number of traces, sorted by trace position, spacing DX. Computation by the original Fortran." );
  pdef->addParam( "apx", "Migration algorithm", NUM_VALUES_FIXED );
  pdef->addValue( "12", VALTYPE_NUMBER, ">= 11: multi-velocity PSPI with APX-10 reference velocities" );
  pdef->addParam( "dx", "Trace spacing within a record (DDD)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "DX" );
  pdef->addParam( "ddz", "Depth sample rate of the image", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "DDZ" );
  pdef->addParam( "nnz", "Depth samples of the image", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NNZ" );
  pdef->addParam( "ltm", "Extra depth samples", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "LTM" );
  pdef->addParam( "itzr", "Depth samples per continuation step", NUM_VALUES_FIXED );
  pdef->addValue( "5", VALTYPE_NUMBER, "ITZR" );
  pdef->addParam( "frq1", "Low frequency of data", NUM_VALUES_FIXED );
  pdef->addValue( "5", VALTYPE_NUMBER, "Hz" );
  pdef->addParam( "frq2", "High frequency of data", NUM_VALUES_FIXED );
  pdef->addValue( "60", VALTYPE_NUMBER, "Hz" );
  pdef->addParam( "oty", "Shift of velocity field (DPN)", NUM_VALUES_FIXED, "Section trace i (1 = first trace of the section) uses velocity DPN i+OTY" );
  pdef->addValue( "0", VALTYPE_NUMBER, "OTY" );
  pdef->addParam( "mute_angle", "Mute angle [deg] (MUTE)", NUM_VALUES_FIXED,
                  "Stack mute: zero above depth GZ + (x/xmax)*x*tan(angle), 20-sample ramp (x = distance to the node). 0 = mute above the node only" );
  pdef->addValue( "30", VALTYPE_NUMBER, "Degrees" );
  pdef->addParam( "lev", "Number of VXT horizons used (LEV1)", NUM_VALUES_FIXED, "0 = all" );
  pdef->addValue( "0", VALTYPE_NUMBER, "LEV1" );
  pdef->addParam( "record_header", "Header that identifies a record", NUM_VALUES_FIXED, "Default: ffid, else spk_rec" );
  pdef->addValue( "ffid", VALTYPE_STRING, "Header name" );
  pdef->addParam( "geometry", "Where the geometry comes from", NUM_VALUES_FIXED );
  pdef->addValue( "headers", VALTYPE_OPTION );
  pdef->addOption( "headers", "Trace positions, source (node) position and depth from trace headers (as ALPH1 = 0 in SEISPAK)" );
  pdef->addOption( "parameters", "ALPH1/GAM1, ALPH2/ALPH3 and ZZY" );
  pdef->addParam( "hdr_trace_x", "Header with the x position of each trace (shot of an OBN receiver gather)", NUM_VALUES_FIXED );
  pdef->addValue( "sou_x", VALTYPE_STRING, "Header name" );
  pdef->addParam( "hdr_trace_y", "Header with the y position of each trace (optional)", NUM_VALUES_FIXED );
  pdef->addValue( "sou_y", VALTYPE_STRING, "Header name" );
  pdef->addParam( "hdr_source_x", "Header with the x position of the source point (node of an OBN receiver gather)", NUM_VALUES_FIXED );
  pdef->addValue( "rec_x", VALTYPE_STRING, "Header name" );
  pdef->addParam( "hdr_source_y", "Header with the y position of the source point (optional)", NUM_VALUES_FIXED );
  pdef->addValue( "rec_y", VALTYPE_STRING, "Header name" );
  pdef->addParam( "hdr_source_z", "Header with the node depth (SEISPAK GZ/GZEL; absolute value used), e.g. rec_wdep", NUM_VALUES_FIXED );
  pdef->addValue( "rec_z", VALTYPE_STRING, "Header name" );
  pdef->addParam( "alph1", "geometry parameters: trace position of the source in record 1 (1 = first trace)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "ALPH1" );
  pdef->addParam( "gam1", "geometry parameters: change of ALPH1 per record", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "GAM1" );
  pdef->addParam( "alph2", "geometry parameters: source depth in record 1", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "ALPH2" );
  pdef->addParam( "alph3", "geometry parameters: change of source depth per record", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "ALPH3" );
  pdef->addParam( "zzy", "geometry parameters: shift of trace 1 per record, in traces", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "ZZY" );
  pdef->addParam( "vxt_file", "File with SEISPAK VXT velocity list (horizons, depth)", NUM_VALUES_FIXED,
                  "DPN 1 = first trace of the section (the leftmost trace position of all records)" );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addParam( "vxt_line", "y-line of the VXT file to use", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "y" );
  pdef->addParam( "vel_file", "SEISPAK file with interval velocity v(z), one trace per DPN (ONL)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addParam( "vel_dz", "Resample the velocity file to this depth step, as SEISPAK onl2dvsvels (DZ)", NUM_VALUES_FIXED,
                  "Interval velocities recomputed preserving the vertical travel time. Default: use the file as it is" );
  pdef->addValue( "", VALTYPE_NUMBER, "DZ" );
  pdef->addParam( "vel_nz", "Number of depth samples after resampling (onl2dvsvels NZ)", NUM_VALUES_FIXED,
                  "Only with vel_dz. Default: as many as the file covers" );
  pdef->addValue( "", VALTYPE_NUMBER, "NZ" );
  pdef->addParam( "focus_x", "Coherence (velocity focussing) analysis of xvspmig1: DPN positions FOR LOR RNI", NUM_VALUES_VARIABLE,
                  "Without focus_x no analysis. At each position the continued fields (receiver x conjugate source) are kept in a window "
                  "around t = 0, depth by depth, summed over the records, then squared, smoothed, converted to time (t2d2), skewed and "
                  "mapped to velocity coherence (tau2v2); picks of the coherence maxima are listed in PREFIX_vtps.txt (vtps)" );
  pdef->addValue( "", VALTYPE_NUMBER, "FOR: first DPN" );
  pdef->addValue( "", VALTYPE_NUMBER, "LOR: last DPN" );
  pdef->addValue( "20", VALTYPE_NUMBER, "RNI: increment" );
  pdef->addParam( "focus_window", "Coherence analysis: time window [ms] (WINDOW)", NUM_VALUES_FIXED );
  pdef->addValue( "1000", VALTYPE_NUMBER, "Window [ms]" );
  pdef->addParam( "focus_vel", "Coherence analysis: velocity range (VMIN VMAX)", NUM_VALUES_FIXED );
  pdef->addValue( "5000", VALTYPE_NUMBER, "VMIN" );
  pdef->addValue( "10000", VALTYPE_NUMBER, "VMAX" );
  pdef->addParam( "focus_eps", "Coherence analysis: threshold of the picked maxima (EPS)", NUM_VALUES_FIXED );
  pdef->addValue( ".1", VALTYPE_NUMBER, "EPS" );
  pdef->addParam( "focus_out", "Coherence analysis: prefix of the outputs PREFIX_tau.sgy, PREFIX_vcoh.sgy, PREFIX_vavg.sgy, PREFIX_vtps.txt", NUM_VALUES_FIXED );
  pdef->addValue( "FOCUS", VALTYPE_STRING, "Prefix" );
  pdef->addParam( "output", "What to output", NUM_VALUES_FIXED );
  pdef->addValue( "stack", VALTYPE_OPTION );
  pdef->addOption( "stack", "Stacked depth section of all records (OFFLINE2); headers trcno, cmp, cmp_x" );
  pdef->addOption( "records", "Each migrated record (not muted), with the input headers (OFFLINE1)" );
}

bool start_exec_mod_spk_xvspmig1_( csExecPhaseEnv* env, csLogWriter* writer ) { return true; }
void cleanup_mod_spk_xvspmig1_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( vars->gather != NULL ) { delete vars->gather; vars->gather = NULL; }
  delete vars; vars = NULL;
}
extern "C" void _params_mod_spk_xvspmig1_( csParamDef* pdef ) { params_mod_spk_xvspmig1_( pdef ); }
extern "C" void _init_mod_spk_xvspmig1_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) { init_mod_spk_xvspmig1_( param, env, writer ); }
extern "C" bool _start_exec_mod_spk_xvspmig1_( csExecPhaseEnv* env, csLogWriter* writer ) { return start_exec_mod_spk_xvspmig1_( env, writer ); }
extern "C" void _exec_mod_spk_xvspmig1_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_xvspmig1_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_xvspmig1_( csExecPhaseEnv* env, csLogWriter* writer ) { cleanup_mod_spk_xvspmig1_( env, writer ); }
