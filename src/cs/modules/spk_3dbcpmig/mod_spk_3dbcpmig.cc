/* Copyright (c) 2026 M. A. Gallotti Guimaraes (COPPE-UFRJ).      */
/* Distributed under the OpenSeaSeis license (see LICENSE).       */
/* Module SPK_3DBCPMIG: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

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
 * Module: SPK_3DBCPMIG
 *
 * 3D prestack depth migration of OBN receiver gathers, from SEISPAK program
 * test3dbcpmig (bcpmig, pfile ptest3dbcpmig, edit exx3dshotmig, DBASE from
 * bcdbasefix). Each record (one node gather, traces = shots) is binned on a local
 * NX x NY grid (DX, DY) centred on the node and aligned with the output grid
 * (XLL/YLL/XLR/YLR), nearest shot per bin, empty bins filled from neighbours
 * within EPS, partial NMO to the bin centre; the recorded wavefield is continued
 * down with 3D multi-velocity PSPI and imaged at the travel time from the node
 * (ITYPE 1: up-going, direct; ITYPE 2: down-going, mirror of the node above the
 * sea surface); angle mute; the image of every record is summed into the output
 * cube NXOUT x NYOUT (OFFLINE2 depth, or OFFLINE3 time).
 * Velocity: SEISPAK ONL v(z) file, one trace per output DPN, one record per output
 * line (VXT -3, ONLVELS). Computation by the original Fortran.
 */

extern "C" {
  void spkbset_( int* nx, int* ny, int* nt, float* apx, int* nnz, int* itzr, float* ddz, float* dt, float* frq1, float* frq2,
                 float* amute, int* nxp, int* nyp, int* nt1, int* nlow, int* nhigh, int* nz, float* zmut, float* hz );
  void spkbv_( float* vdat, int* ntv, int* nxv, int* nyv, float* dzv, int* nxo, int* nyo, int* nz, int* itzr, float* ddz,
               float* dt, int* nt, float* beta, float* vels, float* vrms, int* maxv, float* e, float* x, float* xx );
  void spkbr_( float* d, int* ntd, int* ntr, float* sxh, float* syh, int* igx, int* igy, int* igz, int* ip, float* rp,
               float* vels, float* vrms, float* zmut, float* hz, float* cube, float* wk, float* data, float* tr1, float* yt,
               float* xndx, float* xndx1, float* rr, int* info, int* ifoc, float* tdat );
  void spkbfoc_( float* tdat, float* vels, int* nxout, int* nyout, int* nz, float* dz, int* nt, float* dt, int* jnt,
                 int* ngsnx, int* ngsny, int* iforx, int* ifory, int* irnix, int* irniy, float* vmn, float* vmx,
                 float* otau, float* ovc, float* ovavg, float* zd, float* td, float* s );
  void spkbthr_( int* n );
  void spkbt_( float* x, float* u, float* e, float* s, float* vdat, int* ntv, int* nxv, int* nyv, float* dzv, int* nt,
               float* sr, int* nnz, float* ddz, float* dpnx, float* dpny, float* v );
}

namespace mod_spk_3dbcpmig {
  struct VariableStruct {
    int   itype, nx, ny, nnz, itzr, nt, ifr, ilr, nxout, nyout;
    float dx, dy, ddz, apx, frq1, frq2, amute, gam1, beta, eps;
    float xll, yll, xlr, ylr;
    float sr; int ntIn;
    bool  outTime;
    bool  nodeOv; double nodeX, nodeY, nodeZ;
    std::string recHeader, hdrSx, hdrSy, hdrGx, hdrGy, hdrGz, indataFile;
    int   hmtv, hmrv, ntv; float dzv;
    std::vector<float> vdat;
    csTraceGather* gather;
    // velocity focussing analysis (test3dbcpmig FORX LORX RNIX FORY LORY RNIY WINDOW VMIN VMAX)
    bool  focus; int forx, lorx, rnix, fory, lory, rniy; float window, vmin, vmax;
    std::string tdatFile, focusOut;
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

  // SEASEIS: decimate the velocity grid (extraction to the migration grid): DPN 1, 1+SX, 1+2SX ...; lines 1, 1+SY ...
  void decimateVel( std::vector<float>* vd, int* hmt, int* hmr, int ntv, int sx, int sy ) {
    if( sx <= 1 && sy <= 1 ) return;
    int nx2 = (*hmt - 1)/sx + 1, ny2 = (*hmr - 1)/sy + 1;
    std::vector<float> out( (size_t)nx2*ny2*ntv );
    for( int iy = 0; iy < ny2; iy++ ) for( int ix = 0; ix < nx2; ix++ )
      memcpy( &out[((size_t)iy*nx2 + ix)*ntv], &(*vd)[((size_t)(iy*sy)*(*hmt) + ix*sx)*ntv], ntv*sizeof(float) );
    vd->swap( out ); *hmt = nx2; *hmr = ny2;
  }
  double hval( csTrace* trc, int id ) { return trc->getTraceHeader()->doubleValue( id ); }

  // SEASEIS: small SEG-Y writer (rev 1, IEEE big-endian) for the focussing analysis outputs
  static void putBI( unsigned char* p, int v ) { p[0]=(v>>24)&255; p[1]=(v>>16)&255; p[2]=(v>>8)&255; p[3]=v&255; }
  static void putBS( unsigned char* p, int v ) { p[0]=(v>>8)&255; p[1]=v&255; }
  struct SegyOut {
    FILE* fp; int ns; int count;
    SegyOut() : fp(NULL), ns(0), count(0) {}
    bool open( std::string const& fn, int nsamp, int sampleUs, char const* text ) {
      fp = fopen( fn.c_str(), "wb" ); ns = nsamp; count = 0;
      if( fp == NULL ) return false;
      unsigned char ebc[3200]; memset( ebc, 0x40, 3200 );
      static unsigned char const a2e[128] = {
        0x00,0x01,0x02,0x03,0x37,0x2D,0x2E,0x2F,0x16,0x05,0x25,0x0B,0x0C,0x0D,0x0E,0x0F,0x10,0x11,0x12,0x13,0x3C,0x3D,0x32,0x26,
        0x18,0x19,0x3F,0x27,0x1C,0x1D,0x1E,0x1F,0x40,0x5A,0x7F,0x7B,0x5B,0x6C,0x50,0x7D,0x4D,0x5D,0x5C,0x4E,0x6B,0x60,0x4B,0x61,
        0xF0,0xF1,0xF2,0xF3,0xF4,0xF5,0xF6,0xF7,0xF8,0xF9,0x7A,0x5E,0x4C,0x7E,0x6E,0x6F,0x7C,0xC1,0xC2,0xC3,0xC4,0xC5,0xC6,0xC7,
        0xC8,0xC9,0xD1,0xD2,0xD3,0xD4,0xD5,0xD6,0xD7,0xD8,0xD9,0xE2,0xE3,0xE4,0xE5,0xE6,0xE7,0xE8,0xE9,0xAD,0xE0,0xBD,0x5F,0x6D,
        0x79,0x81,0x82,0x83,0x84,0x85,0x86,0x87,0x88,0x89,0x91,0x92,0x93,0x94,0x95,0x96,0x97,0x98,0x99,0xA2,0xA3,0xA4,0xA5,0xA6,
        0xA7,0xA8,0xA9,0xC0,0x4F,0xD0,0xA1,0x07 };
      for( int i = 0; text[i] != 0 && i < 3200; i++ ) ebc[i] = a2e[ (unsigned char)text[i] & 127 ];
      fwrite( ebc, 1, 3200, fp );
      unsigned char bh[400]; memset( bh, 0, 400 );
      putBI( bh, 1 ); putBS( bh+16, sampleUs ); putBS( bh+20, nsamp ); putBS( bh+24, 5 ); putBS( bh+54, 1 );
      putBS( bh+300, 0x0100 ); putBS( bh+302, 1 );
      fwrite( bh, 1, 400, fp );
      return true;
    }
    void trace( float const* d, int sampleUs, int fldr, int trc, int key, int iline, int xline, double x, double y ) {
      unsigned char th[240]; memset( th, 0, 240 );
      count++;
      putBI( th, count ); putBI( th+4, count ); putBI( th+8, fldr ); putBI( th+12, trc ); putBI( th+20, fldr );
      putBS( th+28, 1 ); putBI( th+36, key ); putBS( th+70, -100 );
      putBS( th+114, ns ); putBS( th+116, sampleUs );
      putBI( th+180, (int)lrint(x*100.0) ); putBI( th+184, (int)lrint(y*100.0) ); putBI( th+188, iline ); putBI( th+192, xline );
      fwrite( th, 1, 240, fp );
      std::vector<unsigned char> b( 4*ns );
      for( int i = 0; i < ns; i++ ) { unsigned int u; memcpy( &u, &d[i], 4 ); putBI( &b[4*i], (int)u ); }
      fwrite( &b[0], 1, b.size(), fp );
    }
    void close() { if( fp != NULL ) fclose( fp ); fp = NULL; }
  };
}
using namespace mod_spk_3dbcpmig;

//*************************************************************************************************
void init_mod_spk_3dbcpmig_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* v = new VariableStruct();
  edef->setVariables( v );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  // defaults of exx3dshotmig
  v->itype = 1; v->nx = 0; v->ny = 0; v->nnz = 0; v->itzr = 5; v->nt = 0; v->ifr = 1; v->ilr = 100000;
  v->nxout = 0; v->nyout = 0; v->dx = 0; v->dy = 0; v->ddz = 0; v->apx = 0; v->frq1 = 0; v->frq2 = 0;
  v->amute = 30; v->gam1 = 0; v->beta = 0; v->eps = 0; v->xll = 0; v->yll = 0; v->xlr = 100; v->ylr = 0;
  v->outTime = false; v->gather = NULL;
  v->focus = false; v->forx = v->lorx = v->rnix = v->fory = v->lory = v->rniy = 0;
  v->window = 1000.0f; v->vmin = 0.0f; v->vmax = 0.0f; v->focusOut = "FOCUS";
  v->hdrSx = "sou_x"; v->hdrSy = "sou_y"; v->hdrGx = "rec_x"; v->hdrGy = "rec_y"; v->hdrGz = "rec_wdep";
  std::string text;
  if( param->exists("itype") ) param->getInt( "itype", &v->itype );
  if( param->exists("nx") )    param->getInt( "nx", &v->nx );
  if( param->exists("ny") )    param->getInt( "ny", &v->ny );
  if( param->exists("dx") )    param->getFloat( "dx", &v->dx );
  if( param->exists("dy") )    param->getFloat( "dy", &v->dy );
  if( param->exists("nnz") )   param->getInt( "nnz", &v->nnz );
  if( param->exists("ddz") )   param->getFloat( "ddz", &v->ddz );
  if( param->exists("itzr") )  param->getInt( "itzr", &v->itzr );
  if( param->exists("nt") )    param->getInt( "nt", &v->nt );
  if( param->exists("apx") )   param->getFloat( "apx", &v->apx );
  if( param->exists("frq1") )  param->getFloat( "frq1", &v->frq1 );
  if( param->exists("frq2") )  param->getFloat( "frq2", &v->frq2 );
  if( param->exists("mute") )  param->getFloat( "mute", &v->amute );
  if( param->exists("elv0") )  param->getFloat( "elv0", &v->gam1 );
  if( param->exists("oty") )   param->getFloat( "oty", &v->beta );
  if( param->exists("eps") )   param->getFloat( "eps", &v->eps );
  if( param->exists("ifr") )   param->getInt( "ifr", &v->ifr );
  if( param->exists("ilr") )   param->getInt( "ilr", &v->ilr );
  if( param->exists("xll") )   param->getFloat( "xll", &v->xll );
  if( param->exists("yll") )   param->getFloat( "yll", &v->yll );
  if( param->exists("xlr") )   param->getFloat( "xlr", &v->xlr );
  if( param->exists("ylr") )   param->getFloat( "ylr", &v->ylr );
  if( param->exists("nxout") ) param->getInt( "nxout", &v->nxout );
  if( param->exists("nyout") ) param->getInt( "nyout", &v->nyout );
  if( param->exists("hdr_shot_x") ) param->getString( "hdr_shot_x", &v->hdrSx );
  if( param->exists("hdr_shot_y") ) param->getString( "hdr_shot_y", &v->hdrSy );
  if( param->exists("hdr_node_x") ) param->getString( "hdr_node_x", &v->hdrGx );
  if( param->exists("hdr_node_y") ) param->getString( "hdr_node_y", &v->hdrGy );
  if( param->exists("hdr_node_z") ) param->getString( "hdr_node_z", &v->hdrGz );
  if( param->exists("indata") ) param->getString( "indata", &v->indataFile );
  v->nodeOv = param->exists("node_xyz");
  if( v->nodeOv ) {
    param->getDouble( "node_xyz", &v->nodeX, 0 );
    param->getDouble( "node_xyz", &v->nodeY, 1 );
    param->getDouble( "node_xyz", &v->nodeZ, 2 );
  }
  if( param->exists("focus_x") || param->exists("focus_y") ) {
    if( !param->exists("focus_x") || !param->exists("focus_y") ) writer->error("Velocity analysis: give both 'focus_x' and 'focus_y'");
    v->focus = true;
    param->getInt( "focus_x", &v->forx, 0 ); param->getInt( "focus_x", &v->lorx, 1 );
    if( param->getNumValues("focus_x") > 2 ) param->getInt( "focus_x", &v->rnix, 2 ); else v->rnix = 1;
    param->getInt( "focus_y", &v->fory, 0 ); param->getInt( "focus_y", &v->lory, 1 );
    if( param->getNumValues("focus_y") > 2 ) param->getInt( "focus_y", &v->rniy, 2 ); else v->rniy = 1;
    if( v->rnix < 1 || v->rniy < 1 || v->lorx < v->forx || v->lory < v->fory || v->forx < 1 || v->fory < 1 )
      writer->error("Velocity analysis: 'focus_x'/'focus_y' must be FIRST LAST STEP, with FIRST >= 1, LAST >= FIRST, STEP >= 1");
    if( param->exists("focus_window") ) param->getFloat( "focus_window", &v->window );
    if( !param->exists("focus_vel") ) writer->error("Velocity analysis: 'focus_vel' (VMIN VMAX) is required");
    param->getFloat( "focus_vel", &v->vmin, 0 ); param->getFloat( "focus_vel", &v->vmax, 1 );
    if( v->vmax <= v->vmin || v->vmin <= 0.0f ) writer->error("Velocity analysis: 'focus_vel' must be VMIN VMAX with 0 < VMIN < VMAX");
    if( param->exists("focus_tdat") ) param->getString( "focus_tdat", &v->tdatFile );
    if( param->exists("focus_out") ) param->getString( "focus_out", &v->focusOut );
  }
  if( param->exists("output") ) {
    param->getString( "output", &text );
    if( !text.compare("time") ) v->outTime = true;
    else if( text.compare("depth") ) writer->error("Unknown option for 'output': %s", text.c_str());
  }
  if( v->itype != 1 && v->itype != 2 ) writer->error("ITYPE %d not supported: 1 = up-going (direct), 2 = down-going (mirror)", v->itype);
  if( v->apx < 11.0f ) writer->error("APX = %g not supported: multi-velocity PSPI only, APX >= 11", v->apx);
  if( v->nx < 4 || v->ny < 1 ) writer->error("NX and NY (local migration grid) are required");
  if( v->dx <= 0.0f ) writer->error("DX is required");
  if( v->dy == 0.0f ) v->dy = v->dx;
  if( v->nnz <= 0 || v->ddz <= 0.0f ) writer->error("NNZ and DDZ are required");
  if( v->itzr < 1 ) writer->error("ITZR must be >= 1");
  if( v->frq2 <= v->frq1 ) writer->error("FRQ2 must be larger than FRQ1");
  if( v->nxout < 1 || v->nyout < 1 ) writer->error("NXOUT and NYOUT (output grid) are required");
  if( v->xll == v->xlr && v->yll == v->ylr ) writer->error("XLL/YLL and XLR/YLR must be different points");

  v->recHeader = "";
  if( param->exists("record_header") ) param->getString( "record_header", &v->recHeader );
  else if( hdef->headerExists("row") ) v->recHeader = "row";
  else if( hdef->headerExists("ffid") ) v->recHeader = "ffid";
  else v->recHeader = "spk_rec";
  std::string const* hn[] = { &v->recHeader, &v->hdrSx, &v->hdrSy, &v->hdrGx, &v->hdrGy, &v->hdrGz };
  for( int i = 0; i < 6; i++ ) if( !hdef->headerExists(*hn[i]) ) writer->error("Trace header '%s' not found", hn[i]->c_str());

  v->ntIn = shdr->numSamples;
  v->sr   = shdr->sampleInt;
  if( v->nt <= 0 || v->nt > v->ntIn ) v->nt = v->ntIn;

  // velocity: SEISPAK ONL (ONLVELS, VXT -3)
  std::string fname;
  if( !param->exists("vel_file") ) writer->error("Parameter 'vel_file' (SEISPAK ONL v(z), as ONLVELS) is required");
  param->getString( "vel_file", &fname );
  std::string err = readSeispak( fname, &v->hmtv, &v->hmrv, &v->ntv, &v->dzv, &v->vdat );
  if( !err.empty() ) { std::string m = "Velocity file '" + fname + "': " + err; writer->error( "%s", m.c_str() ); }
  {
    int sx = 1, sy = 1;
    if( param->exists("vel_step_x") ) param->getInt( "vel_step_x", &sx );
    if( param->exists("vel_step_y") ) param->getInt( "vel_step_y", &sy );
    if( sx < 1 || sy < 1 ) writer->error("vel_step_x / vel_step_y must be >= 1");
    if( sx > 1 || sy > 1 ) {
      int h0 = v->hmtv, r0 = v->hmrv;
      decimateVel( &v->vdat, &v->hmtv, &v->hmrv, v->ntv, sx, sy );
      writer->line("  Velocity grid decimated: every %d DPN and %d line(s): %d x %d -> %d x %d", sx, sy, h0, r0, v->hmtv, v->hmrv);
    }
  }
  writer->line("  Velocity (ONLVELS): %s  %d DPN x %d lines, %d samples, depth step %g", fname.c_str(), v->hmtv, v->hmrv, v->ntv, v->dzv);
  if( v->hmtv < v->nxout ) writer->warning("Velocity grid has %d DPN, fewer than the %d output traces (x): the last DPN is repeated", v->hmtv, v->nxout);
  if( v->hmrv < v->nyout ) writer->warning("Velocity grid has %d lines, fewer than the %d output lines (y): the last line is repeated", v->hmrv, v->nyout);
  {
    int nthr = 0;
    if( param->exists("threads") ) param->getInt( "threads", &nthr );
    spkbthr_( &nthr );
    writer->line("  OpenMP threads (frequency loop): %d", nthr);
  }

  writer->line("  ITYPE %d (%s)  APX %g  ITZR %d  FRQ1/FRQ2 %g %g Hz  MUTE %g deg  ELV0 %g  OTY %g", v->itype,
               v->itype == 1 ? "up-going" : "down-going, mirror", v->apx, v->itzr, v->frq1, v->frq2, v->amute, v->gam1, v->beta);
  writer->line("  Local grid NX %d NY %d  DX %g DY %g  EPS %g;  depth NNZ %d DDZ %g;  NT %d", v->nx, v->ny, v->dx, v->dy, v->eps, v->nnz, v->ddz, v->nt);
  writer->line("  Output grid XLL %.2f YLL %.2f  XLR %.2f YLR %.2f  NXOUT %d NYOUT %d", v->xll, v->yll, v->xlr, v->ylr, v->nxout, v->nyout);
  writer->line("  Records from header '%s' (records %d to %d);  shots '%s' '%s';  node '%s' '%s' depth '%s'", v->recHeader.c_str(),
               v->ifr, v->ilr, v->hdrSx.c_str(), v->hdrSy.c_str(), v->hdrGx.c_str(), v->hdrGy.c_str(), v->hdrGz.c_str());
  if( !v->indataFile.empty() ) writer->line("  Image summed into '%s' (as 'indata')", v->indataFile.c_str());
  if( v->nodeOv ) writer->line("  Node position from parameter node_xyz (DBASE): %.1f %.1f depth %.1f, for every record", v->nodeX, v->nodeY, v->nodeZ);
  if( v->focus ) {
    if( v->lorx > v->nxout || v->lory > v->nyout ) writer->error("Velocity analysis: focus_x/focus_y beyond the output grid (%d x %d)", v->nxout, v->nyout);
    writer->line("  Velocity analysis (focussing): traces %d-%d step %d, lines %d-%d step %d, window %g ms, velocities %g-%g",
                 v->forx, v->lorx, v->rnix, v->fory, v->lory, v->rniy, v->window, v->vmin, v->vmax);
    writer->line("    output '%s_tau.sgy' (OFFLINE4), '%s_vcoh.sgy' (OFFLINE5), '%s_vavg.sgy' (VAVG);  TDAT %s", v->focusOut.c_str(),
                 v->focusOut.c_str(), v->focusOut.c_str(), v->tdatFile.empty() ? "in memory (this run only)" : ("summed into '" + v->tdatFile + "'").c_str());
  }

  v->gather = new csTraceGather( hdef );
  if( !hdef->headerExists(HDR_TRCNO.name) ) hdef->addStandardHeader( HDR_TRCNO.name );
  if( !hdef->headerExists(HDR_ROW.name) )   hdef->addStandardHeader( HDR_ROW.name );
  if( !hdef->headerExists(HDR_COL.name) )   hdef->addStandardHeader( HDR_COL.name );
  if( !hdef->headerExists(HDR_BIN_X.name) ) hdef->addStandardHeader( HDR_BIN_X.name );
  if( !hdef->headerExists(HDR_BIN_Y.name) ) hdef->addStandardHeader( HDR_BIN_Y.name );
  if( v->outTime ) { shdr->numSamples = v->nt; }
  else { shdr->numSamples = v->nnz; shdr->sampleInt = v->ddz; shdr->domain = DOMAIN_XD; }
}

//*************************************************************************************************
void exec_mod_spk_3dbcpmig_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer )
{
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csTraceHeaderDef const* hdef = env->headerDef;
  if( !edef->isLastCall() ) {
    traceGather->moveTracesTo( 0, traceGather->numTraces(), v->gather );
    edef->setTracesAreWaiting();
    return;
  }
  int ntr = v->gather->numTraces();
  if( ntr == 0 ) return;
  int idRec = hdef->headerIndex( v->recHeader );
  std::vector<int> recStart;
  double prev = 0;
  for( int k = 0; k < ntr; k++ ) {
    double h = hval( v->gather->trace(k), idRec );
    if( k == 0 || h != prev ) recStart.push_back( k );
    prev = h;
  }
  int nrecs = (int)recStart.size();
  recStart.push_back( ntr );
  writer->line("  Data: %d records, %d traces, %d samples (%d used)", nrecs, ntr, v->ntIn, v->nt);

  // setup
  int nxp, nyp, nt1, nlow, nhigh, nz;
  float dt = 0.001f*v->sr;
  std::vector<float> zmut( v->nnz + 10 ), hz( v->nnz + 10 );
  spkbset_( &v->nx, &v->ny, &v->nt, &v->apx, &v->nnz, &v->itzr, &v->ddz, &dt, &v->frq1, &v->frq2, &v->amute,
            &nxp, &nyp, &nt1, &nlow, &nhigh, &nz, &zmut[0], &hz[0] );
  if( nhigh < nlow ) writer->error("No frequencies between FRQ1 and FRQ2");
  int nxyout = v->nxout*v->nyout;
  int nxy = nxp*nyp, nws = nhigh - nlow + 1;
  writer->line("  nx ny (FFT) %d %d  nt1 %d  frequency samples %d-%d  depth steps %d (%g m)", nxp, nyp, nt1, nlow, nhigh, nz, v->itzr*v->ddz);

  // velocity
  int maxv = std::max( 4*v->nt, 200000 );
  std::vector<float> vels( (size_t)nxyout*(nz+1) ), vrms( maxv + 10 ), e( 2*v->nt + maxv + 4*nz + 10, 0.0f ), x( v->nnz + nz + 10 ), xx( 2*v->ntv + 10 );
  spkbv_( &v->vdat[0], &v->ntv, &v->hmtv, &v->hmrv, &v->dzv, &v->nxout, &v->nyout, &nz, &v->itzr, &v->ddz, &dt, &v->nt, &v->beta,
          &vels[0], &vrms[0], &maxv, &e[0], &x[0], &xx[0] );
  writer->line("  Mean rms velocity (used for the partial NMO in the binning): %.0f m/s at t = 0, %.0f m/s at t = %g ms", vrms[0], vrms[v->nt-1], (v->nt-1)*v->sr);

  // output cube (indata), nxyout x nz*itzr
  int nzz = nz*v->itzr;
  std::vector<float> cube( (size_t)nxyout*nzz, 0.0f );
  if( !v->indataFile.empty() ) {
    FILE* fp = fopen( v->indataFile.c_str(), "rb" );
    if( fp != NULL ) {
      size_t n = fread( &cube[0], sizeof(float), cube.size(), fp );
      fclose( fp );
      if( n != cube.size() ) writer->error("'%s' has %d values, expected %d (NXOUT*NYOUT*NZ*ITZR): wrong grid?", v->indataFile.c_str(), (int)n, (int)cube.size());
      writer->line("  Image read from '%s' (summing)", v->indataFile.c_str());
    }
    else writer->line("  '%s' does not exist: new image", v->indataFile.c_str());
  }

  // velocity focussing analysis: TDAT(jnt*ngsns, nz), summed node by node like indata
  int ifoc[8] = { 0, 0, 1, 1, 1, 1, 1, 0 };
  int jnt = 1, ngsnx = 0, ngsny = 0;
  if( v->focus ) {
    ngsnx = (v->lorx - v->forx)/v->rnix + 1; ngsny = (v->lory - v->fory)/v->rniy + 1;
    jnt = (int)(v->window/v->sr) + 1;
    if( jnt < 5 ) writer->error("Velocity analysis: window %g ms too short for sample interval %g ms", v->window, v->sr);
    ifoc[0] = ngsnx; ifoc[1] = ngsny; ifoc[2] = v->forx; ifoc[3] = v->fory; ifoc[4] = v->rnix; ifoc[5] = v->rniy; ifoc[6] = jnt;
  }
  int ngsns = ngsnx*ngsny;
  std::vector<float> tdat( std::max( (size_t)1, (size_t)jnt*ngsns*nz ), 0.0f );
  if( v->focus && !v->tdatFile.empty() ) {
    FILE* fp = fopen( v->tdatFile.c_str(), "rb" );
    if( fp != NULL ) {
      size_t n = fread( &tdat[0], sizeof(float), tdat.size(), fp );
      fclose( fp );
      if( n != tdat.size() ) writer->error("'%s' has %d values, expected %d (JNT*points*NZ): different analysis parameters?", v->tdatFile.c_str(), (int)n, (int)tdat.size());
      writer->line("  TDAT read from '%s' (summing)", v->tdatFile.c_str());
    }
    else writer->line("  '%s' does not exist: new TDAT", v->tdatFile.c_str());
  }
  if( v->focus ) writer->line("  Velocity analysis: %d x %d points, JNT %d tau samples (%g ms)", ngsnx, ngsny, jnt, (jnt-1)*v->sr);

  // work space (as bcpmig: t arrays, y/z planes, a b c d, real arrays)
  size_t nxyws = (size_t)nxy*nws;
  size_t id = 1 + 2*(size_t)nxy + 2*nxyws + 3*(size_t)nxy;
  size_t nwk = 2*(id + nxy) + 1 + 4*(size_t)nxy + 2*(size_t)std::max(nxy, nxyout) + 20000;
  std::vector<float> wk( nwk, 0.0f ), data( (size_t)nxy*v->itzr + 10 ), tr1( 4*std::max(v->nt, v->ntIn) + 40 ),
                     yt( 20*nt1 + 1000 ), xndx( nxy + 10 ), xndx1( nxy + 10 ), rr( nxy + 10 );
  writer->line("  Memory %.0f MB", 4.0*(nwk + cube.size() + vels.size())/1.0e6);
  int idSx = hdef->headerIndex(v->hdrSx), idSy = hdef->headerIndex(v->hdrSy);
  int idGx = hdef->headerIndex(v->hdrGx), idGy = hdef->headerIndex(v->hdrGy), idGz = hdef->headerIndex(v->hdrGz);
  int ip[12] = { nxp, nyp, v->nt, nt1, nlow, nhigh, nz, v->itzr, v->nxout, v->nyout, v->itype, 0 };
  float rp[12] = { v->dx, v->dy, dt, v->ddz, v->apx, v->gam1, v->xll, v->yll, v->xlr, v->ylr, v->eps, 0 };
  int ntd = v->ntIn, nmig = 0;
  for( int irec = 1; irec <= nrecs; irec++ ) {
    if( irec < std::max(1, v->ifr) || irec > v->ilr ) continue;
    int k0 = recStart[irec-1], nr = recStart[irec] - k0;
    std::vector<float> d( (size_t)ntd*nr ), sx( nr ), sy( nr );
    for( int k = 0; k < nr; k++ ) {
      csTrace* t = v->gather->trace(k0+k);
      memcpy( &d[(size_t)k*ntd], t->getTraceSamples(), ntd*sizeof(float) );
      sx[k] = (float)hval(t,idSx); sy[k] = (float)hval(t,idSy);
    }
    // DBASE (bcdbasefix): GX GY GZEL of the first trace of the record
    csTrace* t1 = v->gather->trace(k0);
    int igx = (int)lrint( hval(t1,idGx) ), igy = (int)lrint( hval(t1,idGy) ), igz = (int)lrint( hval(t1,idGz) );
    if( v->nodeOv ) { igx = (int)lrint(v->nodeX); igy = (int)lrint(v->nodeY); igz = (int)lrint(v->nodeZ); }   // DBASE given
    int info[10] = {0,0,0,0,0,0,0,0,0,0};
    for( size_t j = 0; j < wk.size(); j++ ) wk[j] = 0.0f;
    spkbr_( &d[0], &ntd, &nr, &sx[0], &sy[0], &igx, &igy, &igz, ip, rp, &vels[0], &vrms[0], &zmut[0], &hz[0], &cube[0],
            &wk[0], &data[0], &tr1[0], &yt[0], &xndx[0], &xndx1[0], &rr[0], info, ifoc, &tdat[0] );
    if( info[2] == 0 ) {
      writer->warning("Record %d (node %d %d, depth %d): node outside the local grid, record skipped", irec, igx, igy, igz);
      continue;
    }
    nmig++;
    writer->line("  Record %d: %d traces, node %d %d depth %d, %d live bins of %d, local grid at output trace %d line %d",
                 irec, nr, igx, igy, igz, info[3], nxy, info[0], info[1]);
  }
  writer->line("  %d records migrated", nmig);
  if( !v->indataFile.empty() ) {
    FILE* fp = fopen( v->indataFile.c_str(), "wb" );
    if( fp == NULL ) writer->error("Cannot write '%s'", v->indataFile.c_str());
    fwrite( &cube[0], sizeof(float), cube.size(), fp );
    fclose( fp );
    writer->line("  Image written to '%s'", v->indataFile.c_str());
  }
  if( v->focus ) {
    if( !v->tdatFile.empty() ) {
      FILE* fp = fopen( v->tdatFile.c_str(), "wb" );
      if( fp == NULL ) writer->error("Cannot write '%s'", v->tdatFile.c_str());
      fwrite( &tdat[0], sizeof(float), tdat.size(), fp );
      fclose( fp );
      writer->line("  TDAT written to '%s'", v->tdatFile.c_str());
    }
    // xvelcoh for every analysis point -> OFFLINE4 (tau deviation), OFFLINE5 (velocity coherence), VAVG
    int nt = v->nt, ntz = std::max( nt, nz );
    float dzf = v->itzr*v->ddz;
    std::vector<float> otau( (size_t)nt*jnt*ngsns ), ovc( (size_t)nt*jnt*ngsns ), ovavg( (size_t)nt*ngsns );
    std::vector<float> zd( (size_t)nz*jnt + 10 ), td( (size_t)(ntz + nt)*jnt + 10 ), ss( 6*(size_t)(nt + nz) + 1000 );
    spkbfoc_( &tdat[0], &vels[0], &v->nxout, &v->nyout, &nz, &dzf, &nt, &dt, &jnt, &ngsnx, &ngsny, &v->forx, &v->fory,
              &v->rnix, &v->rniy, &v->vmin, &v->vmax, &otau[0], &ovc[0], &ovavg[0], &zd[0], &td[0], &ss[0] );
    int iisr = std::max( 1, (int)(8.0f/v->sr + 0.1f) );
    int mjnt = jnt/2 + 1, sus = (int)lrint( v->sr*1000.0 );
    float dv = (v->vmax - v->vmin)/(jnt - 1);
    double rr0 = sqrt( (v->xlr-v->xll)*(v->xlr-v->xll) + (v->ylr-v->yll)*(v->ylr-v->yll) );
    double ux0 = (v->xlr-v->xll)/rr0, uy0 = (v->ylr-v->yll)/rr0;
    SegyOut ft, fv, fa;
    std::string base = v->focusOut;
    if( !ft.open( base + "_tau.sgy", nt, sus, "C 1 SPK_3DBCPMIG VELOCITY FOCUSSING - TAU DEVIATION (OFFLINE4)  FFID = POINT, OFFSET = TAU [MS]" ) ||
        !fv.open( base + "_vcoh.sgy", nt, sus, "C 1 SPK_3DBCPMIG VELOCITY FOCUSSING - VELOCITY COHERENCE (OFFLINE5)  FFID = POINT, OFFSET = VELOCITY" ) ||
        !fa.open( base + "_vavg.sgy", nt, sus, "C 1 SPK_3DBCPMIG VELOCITY FOCUSSING - AVERAGE VELOCITY (VAVG), ONE TRACE PER POINT" ) )
      writer->error("Cannot write the velocity analysis files '%s_*.sgy'", base.c_str());
    for( int iy = 1; iy <= ngsny; iy++ ) {
      for( int ix = 1; ix <= ngsnx; ix++ ) {
        int ii = (iy-1)*ngsnx + ix;
        int dpnx = (ix-1)*v->rnix + v->forx, dpny = (iy-1)*v->rniy + v->fory;
        double bx = v->xll + (dpnx-1)*ux0*v->dx - (dpny-1)*uy0*v->dy, by = v->yll + (dpnx-1)*uy0*v->dx + (dpny-1)*ux0*v->dy;
        for( int j = 1; j <= jnt; j++ ) {
          if( (j-1) % iisr != 0 ) continue;
          int jtrc = (j-1)/iisr + 1;
          size_t o = (size_t)(ii-1)*nt*jnt + (size_t)(j-1)*nt;
          ft.trace( &otau[o], sus, ii, jtrc, (int)lrint((j-mjnt)*v->sr), dpny, dpnx, bx, by );
          fv.trace( &ovc[o], sus, ii, jtrc, (int)lrint(v->vmin + (j-1)*dv), dpny, dpnx, bx, by );
        }
        fa.trace( &ovavg[(size_t)(ii-1)*nt], sus, ii, 1, 0, dpny, dpnx, bx, by );
      }
    }
    ft.close(); fv.close(); fa.close();
    writer->line("  Velocity analysis written: %d points x %d traces (every %d tau samples) to '%s_tau.sgy' and '%s_vcoh.sgy', VAVG to '%s_vavg.sgy'",
                 ngsns, (jnt-1)/iisr + 1, iisr, base.c_str(), base.c_str(), base.c_str());
  }

  // output: OFFLINE2 (depth) or OFFLINE3 (time); record = line (y), trace = x
  int nout = v->outTime ? v->nt : v->nnz;
  traceGather->createTraces( 0, nxyout, hdef, nout );
  double r = sqrt( (v->xlr-v->xll)*(v->xlr-v->xll) + (v->ylr-v->yll)*(v->ylr-v->yll) );
  double ux = (v->xlr-v->xll)/r, uy = (v->ylr-v->yll)/r;
  int idTrc = hdef->headerIndex(HDR_TRCNO.name), idRow = hdef->headerIndex(HDR_ROW.name), idCol = hdef->headerIndex(HDR_COL.name);
  int idBx = hdef->headerIndex(HDR_BIN_X.name), idBy = hdef->headerIndex(HDR_BIN_Y.name);
  std::vector<float> xz( v->nnz + 10 ), u( v->nt + 10 ), ee( v->nt + 10 ), ss( v->nnz + 10 ), vv( 2*v->ntv + 10 );
  for( int k = 1; k <= v->nyout; k++ ) {
    for( int i = 1; i <= v->nxout; i++ ) {
      int irec = (k-1)*v->nxout + i;
      csTrace* trc = traceGather->trace( irec-1 );
      csTraceHeader* th = trc->getTraceHeader();
      th->setIntValue( idTrc, irec );
      th->setIntValue( idRow, k );
      th->setIntValue( idCol, i );
      th->setDoubleValue( idBx, v->xll + (i-1)*ux*v->dx - (k-1)*uy*v->dy );
      th->setDoubleValue( idBy, v->yll + (i-1)*uy*v->dx + (k-1)*ux*v->dy );
      float* out = trc->getTraceSamples();
      for( int j = 0; j < v->nnz; j++ ) xz[j] = ( j < nzz ) ? cube[(size_t)j*nxyout + irec-1] : 0.0f;
      if( !v->outTime ) memcpy( out, &xz[0], v->nnz*sizeof(float) );
      else {
        float dpnx = (float)i, dpny = (float)k;
        spkbt_( &xz[0], &u[0], &ee[0], &ss[0], &v->vdat[0], &v->ntv, &v->hmtv, &v->hmrv, &v->dzv, &v->nt, &v->sr, &v->nnz, &v->ddz,
                &dpnx, &dpny, &vv[0] );
        memcpy( out, &u[0], v->nt*sizeof(float) );
      }
    }
  }
  writer->line("  Output: %d lines x %d traces, %s", v->nyout, v->nxout, v->outTime ? "time (OFFLINE3)" : "depth (OFFLINE2)");
  v->gather->freeAllTraces();
}

//********************************************************************************
void params_mod_spk_3dbcpmig_( csParamDef* pdef ) {
  pdef->setModule( "SPK_3DBCPMIG", "3D prestack depth migration of OBN receiver gathers (SEISPAK test3dbcpmig)",
                   "Each record (node gather) is binned on a local NX x NY grid centred on the node and aligned with the output grid, "
                   "continued down with 3D multi-velocity PSPI and imaged at the travel time from the node (ITYPE 1 up-going, ITYPE 2 "
                   "down-going with the mirror of the node). The images of all records are summed into the NXOUT x NYOUT output cube. "
                   "Computation by the original Fortran." );
  pdef->addParam( "itype", "1 = up-going (direct), 2 = down-going (mirror)", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "ITYPE" );
  pdef->addParam( "apx", "Migration algorithm (>= 11: multi-velocity PSPI, APX-10 velocities)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "APX" );
  pdef->addParam( "nx", "Local grid: inline (x) bins", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NX" );
  pdef->addParam( "ny", "Local grid: crossline (y) bins", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NY" );
  pdef->addParam( "dx", "Bin size in x (GGD inline)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "DX" );
  pdef->addParam( "dy", "Bin size in y (default DX)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "DY" );
  pdef->addParam( "eps", "Maximum distance for filling an empty bin from its neighbours", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "EPS" );
  pdef->addParam( "nnz", "Depth samples", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NNZ" );
  pdef->addParam( "ddz", "Depth sample rate", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "DDZ" );
  pdef->addParam( "itzr", "Depth samples per continuation step", NUM_VALUES_FIXED );
  pdef->addValue( "5", VALTYPE_NUMBER, "ITZR" );
  pdef->addParam( "nt", "Time samples used (default: all)", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "NT" );
  pdef->addParam( "frq1", "Low frequency [Hz]", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "FRQ1" );
  pdef->addParam( "frq2", "High frequency [Hz]", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "FRQ2" );
  pdef->addParam( "mute", "Mute angle [deg]", NUM_VALUES_FIXED );
  pdef->addValue( "30", VALTYPE_NUMBER, "MUTE" );
  pdef->addParam( "elv0", "ELV0: depth where the continuation with the velocity field starts (kstart = ELV0/dz+1)", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "ELV0" );
  pdef->addParam( "oty", "Shift of the velocity DPN: output trace i uses DPN i+OTY", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "OTY" );
  pdef->addParam( "ifr", "First record to migrate", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "IFR" );
  pdef->addParam( "ilr", "Last record to migrate", NUM_VALUES_FIXED );
  pdef->addValue( "100000", VALTYPE_NUMBER, "ILR" );
  pdef->addParam( "xll", "Output grid: x of the lower left corner (origin)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "XLL" );
  pdef->addParam( "yll", "Output grid: y of the lower left corner (origin)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "YLL" );
  pdef->addParam( "xlr", "Output grid: x of a point along the inline (x) axis", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "XLR" );
  pdef->addParam( "ylr", "Output grid: y of a point along the inline (x) axis", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "YLR" );
  pdef->addParam( "nxout", "Output grid: traces in x", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NXOUT" );
  pdef->addParam( "nyout", "Output grid: lines in y", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NYOUT" );
  pdef->addParam( "vel_file", "SEISPAK ONL v(z): one trace per output DPN, one record per output line (ONLVELS, VXT -3)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addParam( "vel_step_x", "Velocity grid: DPN per output trace (extraction: output trace i uses DPN (i-1)*STEP+1)", NUM_VALUES_FIXED,
                  "Ex.: model at 25 m, migration grid at 50 m, same origin and azimuth: 2" );
  pdef->addValue( "1", VALTYPE_NUMBER, "Step" );
  pdef->addParam( "vel_step_y", "Velocity grid: lines per output line", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "Step" );
  pdef->addParam( "record_header", "Header that identifies a record (node gather)", NUM_VALUES_FIXED, "Default: row, else ffid, else spk_rec" );
  pdef->addValue( "row", VALTYPE_STRING, "Header name" );
  pdef->addParam( "hdr_shot_x", "Shot x (SEISPAK SX)", NUM_VALUES_FIXED );
  pdef->addValue( "sou_x", VALTYPE_STRING, "Header name" );
  pdef->addParam( "hdr_shot_y", "Shot y (SEISPAK SY)", NUM_VALUES_FIXED );
  pdef->addValue( "sou_y", VALTYPE_STRING, "Header name" );
  pdef->addParam( "hdr_node_x", "Node x (SEISPAK GX, first trace of the record)", NUM_VALUES_FIXED );
  pdef->addValue( "rec_x", VALTYPE_STRING, "Header name" );
  pdef->addParam( "hdr_node_y", "Node y (SEISPAK GY)", NUM_VALUES_FIXED );
  pdef->addValue( "rec_y", VALTYPE_STRING, "Header name" );
  pdef->addParam( "hdr_node_z", "Node depth (SEISPAK GZEL; absolute value used)", NUM_VALUES_FIXED );
  pdef->addValue( "rec_wdep", VALTYPE_STRING, "Header name" );
  pdef->addParam( "node_xyz", "Node x y depth for every record, instead of the headers (as editing the DBASE line)", NUM_VALUES_FIXED,
                  "To reproduce a SEISPAK run whose DBASE had other coordinates" );
  pdef->addValue( "", VALTYPE_NUMBER, "x" );
  pdef->addValue( "", VALTYPE_NUMBER, "y" );
  pdef->addValue( "", VALTYPE_NUMBER, "depth" );
  pdef->addParam( "indata", "Binary file with the image being summed (as the SEISPAK 'indata' link): read if it exists, written at the end", NUM_VALUES_FIXED,
                  "Lets several runs (nodes) add into the same cube" );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addParam( "threads", "Number of CPU threads (OpenMP, frequencies migrated in parallel)", NUM_VALUES_FIXED,
                  "Default 0: OMP_NUM_THREADS if set, else all cores. Same result for any number of threads" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Threads" );
  pdef->addParam( "focus_x", "Velocity analysis (focussing, test3dbcpmig FORX LORX RNIX): analysis points along x, in output traces", NUM_VALUES_VARIABLE,
                  "Without focus_x/focus_y no analysis is done. At each point the continued wavefield is kept in a window around the arrival time, "
                  "depth by depth (TDAT, summed over the nodes), then converted to tau-deviation and velocity-coherence panels (xvelcoh)" );
  pdef->addValue( "", VALTYPE_NUMBER, "FORX: first output trace" );
  pdef->addValue( "", VALTYPE_NUMBER, "LORX: last output trace" );
  pdef->addValue( "1", VALTYPE_NUMBER, "RNIX: increment" );
  pdef->addParam( "focus_y", "Velocity analysis: analysis points along y, in output lines (FORY LORY RNIY)", NUM_VALUES_VARIABLE );
  pdef->addValue( "", VALTYPE_NUMBER, "FORY: first output line" );
  pdef->addValue( "", VALTYPE_NUMBER, "LORY: last output line" );
  pdef->addValue( "1", VALTYPE_NUMBER, "RNIY: increment" );
  pdef->addParam( "focus_window", "Velocity analysis: tau window centred on the arrival time [ms] (WINDOW)", NUM_VALUES_FIXED );
  pdef->addValue( "1000", VALTYPE_NUMBER, "Window [ms]" );
  pdef->addParam( "focus_vel", "Velocity analysis: velocity range of the coherence panels (VMIN VMAX)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "VMIN" );
  pdef->addValue( "", VALTYPE_NUMBER, "VMAX" );
  pdef->addParam( "focus_tdat", "Velocity analysis: binary file with TDAT, read if it exists and written at the end (as 'indata')", NUM_VALUES_FIXED,
                  "Lets several runs (nodes) add into the same analysis; delete it to start a new one" );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addParam( "focus_out", "Velocity analysis: prefix of the output SEG-Y files PREFIX_tau.sgy, PREFIX_vcoh.sgy, PREFIX_vavg.sgy", NUM_VALUES_FIXED );
  pdef->addValue( "FOCUS", VALTYPE_STRING, "Prefix" );
  pdef->addParam( "output", "Output cube", NUM_VALUES_FIXED );
  pdef->addValue( "depth", VALTYPE_OPTION );
  pdef->addOption( "depth", "Depth cube (OFFLINE2): NNZ samples at DDZ" );
  pdef->addOption( "time", "Time cube (OFFLINE3): NT samples, converted with the velocity of each DPN" );
}

bool start_exec_mod_spk_3dbcpmig_( csExecPhaseEnv* env, csLogWriter* writer ) { return true; }
void cleanup_mod_spk_3dbcpmig_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( v->gather != NULL ) { delete v->gather; v->gather = NULL; }
  delete v; v = NULL;
}
extern "C" void _params_mod_spk_3dbcpmig_( csParamDef* pdef ) { params_mod_spk_3dbcpmig_( pdef ); }
extern "C" void _init_mod_spk_3dbcpmig_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) { init_mod_spk_3dbcpmig_( param, env, writer ); }
extern "C" bool _start_exec_mod_spk_3dbcpmig_( csExecPhaseEnv* env, csLogWriter* writer ) { return start_exec_mod_spk_3dbcpmig_( env, writer ); }
extern "C" void _exec_mod_spk_3dbcpmig_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_3dbcpmig_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_3dbcpmig_( csExecPhaseEnv* env, csLogWriter* writer ) { cleanup_mod_spk_3dbcpmig_( env, writer ); }
