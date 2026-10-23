/* Copyright (c) 2026 M. A. Gallotti Guimaraes (COPPE-UFRJ).      */
/* Distributed under the OpenSeaSeis license (see LICENSE).       */
/* Module SPK_SHOTBIN: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

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
 * Module: SPK_SHOTBIN
 *
 * SEISPAK program XXSHOTBIN (process BIN): the shots of one receiver gather
 * (OBN/PRM) are placed on a regular, rotated grid of NX x NY bins defined by
 * the lower left corner (XLL,YLL), the direction to (XLR,YLR) and DX, DY.
 * For each bin the trace whose source is closest to the bin centre (and
 * within EPS) is kept; with PMO 1 it is moved out from its own offset to the
 * offset of the bin centre (hyperbolic, v(t) from VEL). Empty bins take the
 * nearest trace of the 8 neighbouring bins if within EPS, otherwise a zero trace.
 * The SPL (spline) and INT3 (velocity interpolation) routines are the original
 * Fortran (fortran/shotbin_sub.f); the main program (jhh I/O, TABLE1) is here.
 *
 * Changes to the original, marked SEASEIS:
 *  - Coordinates in double precision (Fortran REAL: ~0.25 m at y = 3e6 m).
 *  - Neighbour search for empty bins checks the grid limits (the original
 *    wrapped around to the other side of the grid at the x edges and could
 *    read outside the index array at the y edges).
 *  - t0 for the bin-centre offset is computed for each trace with the same
 *    three fixed-point steps (in the original the steps carried over from the
 *    previous bin, so the result depended on the bin order).
 *  - GXBIN/GYBIN are used for the moveout (the original only printed them).
 *  - Dead bins take the headers of the first input trace (the original left the
 *    headers of the previous bin).
 *  - INDEX/XNDX (never written or read in the original) and COORDS/BCOORDS/BINPLOT
 *    are replaced by the optional ASCII QC file (qc_file).
 */

extern "C" {
  void spl_( float* y, float* x, float* f, float* s, int* n, int* m );
  void int3_( float* f, float* x, float* y, float* x1, float* dx, int* n, int* m );
}

namespace mod_spk_shotbin {
  struct VariableStruct {
    double xll, yll, xlr, ylr;
    double dx, dy, eps;
    int    nx, ny;
    int    hmt;              // >1: records of NX traces (NY records); else one record
    int    pmo;
    float  gxbin, gybin;
    std::vector<float> velPairs;
    std::string qcFile;
    int    nt;
    float  sr;
    // header indices
    int idSx, idSy, idGx, idGy, idGz;
    int idRang, idOffset, idRnum, idTnum;
    int idBinX, idBinY, idRow, idCol, idDist, idFlag;
    csTraceGather* gather;
  };

  // First existing header name of a list
  static int findHeader( csTraceHeaderDef* hdef, char const* const* names, int n, std::string* found ) {
    for( int i = 0; i < n; i++ ) {
      if( hdef->headerExists( names[i] ) ) { *found = names[i]; return hdef->headerIndex( names[i] ); }
    }
    return -1;
  }
}
using namespace mod_spk_shotbin;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_spk_shotbin_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* vars = new VariableStruct();
  edef->setVariables( vars );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );
  vars->gather = NULL;

  std::string text;
  char const* req[] = { "xll", "yll", "xlr", "ylr", "dx", "dy", "nx", "ny", "eps" };
  for( int i = 0; i < 9; i++ ) {
    if( !param->exists( req[i] ) ) writer->error("Parameter '%s' is required", req[i]);
  }
  param->getDouble( "xll", &vars->xll );
  param->getDouble( "yll", &vars->yll );
  param->getDouble( "xlr", &vars->xlr );
  param->getDouble( "ylr", &vars->ylr );
  param->getDouble( "dx", &vars->dx );
  param->getDouble( "dy", &vars->dy );
  param->getInt( "nx", &vars->nx );
  param->getInt( "ny", &vars->ny );
  param->getDouble( "eps", &vars->eps );
  vars->hmt = 1;
  if( param->exists("hmt") ) { float h; param->getFloat( "hmt", &h ); vars->hmt = (int)lrintf(h); }
  vars->pmo = 0;
  if( param->exists("pmo") ) { float p; param->getFloat( "pmo", &p ); vars->pmo = (p == 1.0f) ? 1 : 0; }
  vars->gxbin = 0.0f; vars->gybin = 0.0f;
  if( param->exists("gxbin") ) param->getFloat( "gxbin", &vars->gxbin );
  if( param->exists("gybin") ) param->getFloat( "gybin", &vars->gybin );
  if( param->exists("qc_file") ) param->getString( "qc_file", &vars->qcFile );

  if( vars->nx < 1 || vars->ny < 1 ) writer->error("NX and NY must be >= 1");
  if( vars->dx <= 0.0 || vars->dy <= 0.0 ) writer->error("DX and DY must be > 0");
  if( vars->eps <= 0.0 ) writer->error("EPS must be > 0");
  double r = sqrt( (vars->xlr-vars->xll)*(vars->xlr-vars->xll) + (vars->ylr-vars->yll)*(vars->ylr-vars->yll) );
  if( r <= 0.0 ) writer->error("(XLR,YLR) must differ from (XLL,YLL): it gives the direction of the x axis");

  // Velocity-time pairs (all values of all 'vel' lines)
  int nLines = param->getNumLines( "vel" );
  for( int il = 0; il < nLines; il++ ) {
    int nv = param->getNumValues( "vel", il );
    for( int iv = 0; iv < nv; iv++ ) {
      float v; param->getFloatAtLine( "vel", &v, il, iv );
      vars->velPairs.push_back( v );
    }
  }
  if( vars->pmo == 1 ) {
    int nvel = (int)vars->velPairs.size();
    if( nvel < 2 || nvel % 2 != 0 ) writer->error("PMO 1 requires VEL: pairs of velocity and time [ms] (v1 t1 v2 t2 ...)");
    for( int j = 0; j < nvel/2; j++ ) {
      if( vars->velPairs[2*j] <= 0.0f ) writer->error("VEL: velocities must be > 0");
      if( j > 0 && vars->velPairs[2*j+1] <= vars->velPairs[2*j-1] ) writer->error("VEL: times must increase");
    }
  }

  vars->nt = shdr->numSamples;
  vars->sr = shdr->sampleInt;
  if( vars->pmo == 1 && vars->nt < 4 ) writer->error("Spline interpolation (PMO 1) needs at least 4 samples");

  //---------------------------------------------
  // Headers: source x/y (HDRS), receiver x/y/z
  std::string sxName, syName, gxName, gyName, gzName;
  if( param->exists("hdrs") ) {
    param->getString( "hdrs", &sxName, 0 );
    param->getString( "hdrs", &syName, 1 );
    if( !hdef->headerExists(sxName) ) writer->error("HDRS: header '%s' not found", sxName.c_str());
    if( !hdef->headerExists(syName) ) writer->error("HDRS: header '%s' not found", syName.c_str());
    vars->idSx = hdef->headerIndex( sxName );
    vars->idSy = hdef->headerIndex( syName );
  }
  else {
    char const* nsx[] = { "spk_sx", "sou_x" };
    char const* nsy[] = { "spk_sy", "sou_y" };
    vars->idSx = findHeader( hdef, nsx, 2, &sxName );
    vars->idSy = findHeader( hdef, nsy, 2, &syName );
    if( vars->idSx < 0 || vars->idSy < 0 ) writer->error("Source coordinate headers not found (spk_sx/spk_sy or sou_x/sou_y). Use 'hdrs'");
  }
  if( param->exists("hdr_rec") ) {
    param->getString( "hdr_rec", &gxName, 0 );
    param->getString( "hdr_rec", &gyName, 1 );
    param->getString( "hdr_rec", &gzName, 2 );
    if( !hdef->headerExists(gxName) || !hdef->headerExists(gyName) || !hdef->headerExists(gzName) ) {
      writer->error("HDR_REC: header(s) not found: %s %s %s", gxName.c_str(), gyName.c_str(), gzName.c_str());
    }
    vars->idGx = hdef->headerIndex( gxName );
    vars->idGy = hdef->headerIndex( gyName );
    vars->idGz = hdef->headerIndex( gzName );
  }
  else {
    char const* ngx[] = { "spk_gx", "rec_x" };
    char const* ngy[] = { "spk_gy", "rec_y" };
    char const* ngz[] = { "spk_gz", "spk_gzel", "rec_z", "rec_wdep" };   // GZ, else GZEL (as the original)
    vars->idGx = findHeader( hdef, ngx, 2, &gxName );
    vars->idGy = findHeader( hdef, ngy, 2, &gyName );
    vars->idGz = findHeader( hdef, ngz, 4, &gzName );
    if( vars->idGx < 0 || vars->idGy < 0 ) writer->error("Receiver coordinate headers not found (spk_gx/spk_gy or rec_x/rec_y). Use 'hdr_rec'");
  }
  vars->idRang   = hdef->headerExists("spk_rang") ? hdef->headerIndex("spk_rang") : -1;
  vars->idOffset = hdef->headerExists("offset")   ? hdef->headerIndex("offset")   : -1;
  vars->idRnum   = hdef->headerExists("spk_rnum") ? hdef->headerIndex("spk_rnum") : -1;
  vars->idTnum   = hdef->headerExists("spk_tnum") ? hdef->headerIndex("spk_tnum") : -1;
  vars->idBinX = hdef->addStandardHeader( HDR_BIN_X.name );
  vars->idBinY = hdef->addStandardHeader( HDR_BIN_Y.name );
  vars->idRow  = hdef->addStandardHeader( HDR_ROW.name );
  vars->idCol  = hdef->addStandardHeader( HDR_COL.name );
  if( !hdef->headerExists("shotbin_dist") ) hdef->addHeader( TYPE_FLOAT, "shotbin_dist", "SPK_SHOTBIN: distance source - bin centre [m] (-1: dead bin)" );
  if( !hdef->headerExists("shotbin_flag") ) hdef->addHeader( TYPE_INT, "shotbin_flag", "SPK_SHOTBIN: 0 dead bin, 1 trace of the bin, 2 trace from a neighbouring bin" );
  vars->idDist = hdef->headerIndex( "shotbin_dist" );
  vars->idFlag = hdef->headerIndex( "shotbin_flag" );

  vars->gather = new csTraceGather( hdef );

  writer->line("  Grid: NX %d  NY %d  DX %g  DY %g  EPS %g", vars->nx, vars->ny, vars->dx, vars->dy, vars->eps);
  writer->line("  XLL YLL %.2f %.2f   XLR YLR %.2f %.2f", vars->xll, vars->yll, vars->xlr, vars->ylr);
  writer->line("  Source headers (HDRS): %s %s   Receiver headers: %s %s %s", sxName.c_str(), syName.c_str(),
               gxName.c_str(), gyName.c_str(), vars->idGz >= 0 ? gzName.c_str() : "(none)");
  writer->line("  PMO %d%s", vars->pmo, vars->pmo ? " (partial moveout to bin centre)" : " (no moveout)");
  if( vars->pmo ) {
    for( size_t j = 0; j+1 < vars->velPairs.size(); j += 2 ) writer->line("     v %8.1f  at t %8.1f ms", vars->velPairs[j], vars->velPairs[j+1]);
  }
  writer->line("  Output: %d traces, %s", vars->nx*vars->ny,
               vars->hmt > 1 ? "records of NX traces (one record per y-line)" : "one record");
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_spk_shotbin_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csExecPhaseDef* edef = env->execPhaseDef;
  csTraceHeaderDef const* hdef = env->headerDef;

  if( !edef->isLastCall() ) {
    traceGather->moveTracesTo( 0, traceGather->numTraces(), vars->gather );
    edef->setTracesAreWaiting();
    return;
  }
  int ntrIn = vars->gather->numTraces();
  if( ntrIn == 0 ) return;

  int nx = vars->nx, ny = vars->ny, nt = vars->nt;
  int nbin = nx*ny;
  double dt = 0.001 * vars->sr;

  //--- grid axes (as in the original)
  double r  = sqrt( (vars->xlr-vars->xll)*(vars->xlr-vars->xll) + (vars->ylr-vars->yll)*(vars->ylr-vars->yll) );
  double ux = (vars->xlr-vars->xll)/r;
  double uy = (vars->ylr-vars->yll)/r;
  double dxx = ux*vars->dx, dxy = -uy*vars->dy, dyx = uy*vars->dx, dyy = ux*vars->dy;
  double xll = vars->xll, yll = vars->yll;
  writer->line("  dxx dxy dyx dyy  %f %f %f %f", dxx, dxy, dyx, dyy);
  writer->line("  Corners: LL %.1f %.1f  LR %.1f %.1f  UL %.1f %.1f  UR %.1f %.1f", xll, yll,
               xll+(nx-1)*dxx, yll+(nx-1)*dyx, xll+(ny-1)*dxy, yll+(ny-1)*dyy,
               xll+(nx-1)*dxx+(ny-1)*dxy, yll+(nx-1)*dyx+(ny-1)*dyy);

  //--- velocity v(t) at every sample (INT3), stored as 1/(v*dt)^2 (samples)
  std::vector<float> vel( nt+1, 0.0f );   // 1-based
  if( vars->pmo == 1 ) {
    int nvel = (int)vars->velPairs.size();
    int npairs = nvel/2 + 2;
    std::vector<float> xs( npairs ), ys( npairs );
    ys[0] = vars->velPairs[0]; xs[0] = 0.0f;
    for( int j = 1; j <= nvel/2; j++ ) { ys[j] = vars->velPairs[2*j-2]; xs[j] = vars->velPairs[2*j-1]; }
    ys[npairs-1] = ys[npairs-2];
    xs[npairs-1] = nt * vars->sr;
    if( xs[npairs-2] >= xs[npairs-1] ) xs[npairs-1] = xs[npairs-2] + vars->sr;  // SEASEIS: keep times increasing
    float x1 = 0.0f, sr = vars->sr;
    int n = nt, m = npairs;
    int3_( &vel[1], &xs[0], &ys[0], &x1, &sr, &n, &m );
    for( int j = 1; j <= nt; j++ ) {
      if( !(vel[j] > 0.0f) ) vel[j] = ys[npairs-1];
      vel[j] = (float)( 1.0 / ( (double)vel[j]*dt * (double)vel[j]*dt ) );
    }
  }

  //--- reference receiver (first trace) and check that the input is one receiver gather
  csTraceHeader const* h0 = vars->gather->trace(0)->getTraceHeader();
  double gxRef = h0->doubleValue( vars->idGx );
  double gyRef = h0->doubleValue( vars->idGy );
  double gzRef = ( vars->idGz >= 0 ) ? h0->doubleValue( vars->idGz ) : 0.0;
  double maxRecDiff = 0.0;
  for( int k = 1; k < ntrIn; k++ ) {
    csTraceHeader const* th = vars->gather->trace(k)->getTraceHeader();
    double d = sqrt( pow( th->doubleValue(vars->idGx)-gxRef, 2 ) + pow( th->doubleValue(vars->idGy)-gyRef, 2 ) );
    if( d > maxRecDiff ) maxRecDiff = d;
  }
  if( maxRecDiff > vars->eps ) {
    writer->warning("Receiver position varies by up to %.1f m in the input: SPK_SHOTBIN bins ONE receiver gather "
                    "(output receiver = first trace). Select one receiver before this module.", maxRecDiff);
  }
  // GXBIN/GYBIN: receiver position for the moveout given as bin numbers (SEASEIS: used, see header)
  bool recOverride = ( vars->gxbin > 0.0f );
  double gxOv = 0.0, gyOv = 0.0;
  if( recOverride ) {
    gxOv = xll + (vars->gxbin-1.0)*dxx + (vars->gybin-1.0)*dxy;
    gyOv = yll + (vars->gxbin-1.0)*dyx + (vars->gybin-1.0)*dyy;
    writer->line("  Receiver for moveout (GXBIN %g GYBIN %g): gx = %.2f  gy = %.2f", vars->gxbin, vars->gybin, gxOv, gyOv);
  }
  {
    double xg = ( (gxRef-xll)*ux + (gyRef-yll)*uy )/vars->dx + 1.0;
    double yg = ( -(gxRef-xll)*uy + (gyRef-yll)*ux )/vars->dy + 1.0;
    writer->line("  Receiver: gx %.2f gy %.2f gz %.2f  ->  bin x %.1f  y %.1f", gxRef, gyRef, gzRef, xg, yg);
  }

  //--- pass 1: trace of each bin (closest source within EPS)
  std::vector<double> sxv( ntrIn ), syv( ntrIn );
  std::vector<int>    indx( nbin, -1 ), flag( nbin, 0 );
  std::vector<double> rr( nbin, 1.0e7 );
  std::vector<int>    hits( nbin, 0 );       // sources falling in each bin (BINPLOT 'dots')
  int ndup = 0, nOutside = 0;
  for( int k = 0; k < ntrIn; k++ ) {
    csTraceHeader const* th = vars->gather->trace(k)->getTraceHeader();
    double sx = th->doubleValue( vars->idSx );
    double sy = th->doubleValue( vars->idSy );
    sxv[k] = sx; syv[k] = sy;
    double xxsx = ( (sx-xll)*ux + (sy-yll)*uy )/vars->dx + 1.0;
    double yysy = ( (sx-xll)*(-uy) + (sy-yll)*ux )/vars->dy + 1.0;
    int ixbin = (int)lround( xxsx );
    int iybin = (int)lround( yysy );
    if( ixbin < 1 || ixbin > nx || iybin < 1 || iybin > ny ) { nOutside++; continue; }
    int ixy = (iybin-1)*nx + ixbin - 1;
    hits[ixy] += 1;
    double xc = xll + (ixbin-1)*dxx + (iybin-1)*dxy;
    double yc = yll + (ixbin-1)*dyx + (iybin-1)*dyy;
    double rng = sqrt( (xc-sx)*(xc-sx) + (yc-sy)*(yc-sy) );
    if( rng < rr[ixy] && rng < vars->eps ) {
      if( rr[ixy] < 1.0e6 ) ndup++;
      rr[ixy] = rng;
      indx[ixy] = k;
      flag[ixy] = 1;
    }
  }

  //--- pass 2: empty bins take the closest trace of the 8 neighbours (within EPS)
  // As in the original, a bin filled here can serve the following bins (indx is updated in the loop)
  for( int iy = 1; iy <= ny; iy++ ) {
    for( int ix = 1; ix <= nx; ix++ ) {
      int ii = (iy-1)*nx + ix - 1;
      if( indx[ii] >= 0 ) continue;
      double xc = xll + (ix-1)*dxx + (iy-1)*dxy;
      double yc = yll + (ix-1)*dyx + (iy-1)*dyy;
      for( int i1 = -1; i1 <= 1; i1++ ) {
        for( int i2 = -1; i2 <= 1; i2++ ) {
          if( i1 == 0 && i2 == 0 ) continue;
          int jy = iy + i1, jx = ix + i2;
          if( jy < 1 || jy > ny || jx < 1 || jx > nx ) continue;   // SEASEIS: grid limits
          int k = indx[(jy-1)*nx + jx - 1];
          if( k < 0 ) continue;
          double rng = sqrt( (xc-sxv[k])*(xc-sxv[k]) + (yc-syv[k])*(yc-syv[k]) );
          if( rng < rr[ii] && rng < vars->eps ) {
            rr[ii] = rng;
            indx[ii] = k;
            flag[ii] = 2;
          }
        }
      }
    }
  }

  //--- output traces
  traceGather->createTraces( 0, nbin, hdef, nt );
  std::vector<float> s1( 3*nt + 3 );
  std::vector<int>   it0( nt+1 );
  int ihmt1 = ( vars->hmt > 1 ) ? nx : nbin;
  int nlive = 0, ndead = 0, nfill = 0;
  double avgrng = 0.0;
  for( int iy = 1; iy <= ny; iy++ ) {
    for( int ix = 1; ix <= nx; ix++ ) {
      int ii = (iy-1)*nx + ix - 1;
      double xc = xll + (ix-1)*dxx + (iy-1)*dxy;
      double yc = yll + (ix-1)*dyx + (iy-1)*dyy;
      csTrace* trcOut = traceGather->trace( ii );
      csTraceHeader* hdrOut = trcOut->getTraceHeader();
      float* out = trcOut->getTraceSamples();
      int k = indx[ii];
      if( k >= 0 ) {
        nlive++;
        if( flag[ii] == 2 ) nfill++;
        csTrace* trcIn = vars->gather->trace( k );
        hdrOut->copyFrom( trcIn->getTraceHeader() );
        float* in = trcIn->getTraceSamples();
        avgrng += rr[ii];
        if( vars->pmo == 1 ) {
          double gx = recOverride ? gxOv : hdrOut->doubleValue( vars->idGx );
          double gy = recOverride ? gyOv : hdrOut->doubleValue( vars->idGy );
          double sx = sxv[k], sy = syv[k];
          double rng2 = (gx-sx)*(gx-sx) + (gy-sy)*(gy-sy);     // offset^2 of the trace
          double rng1 = (gx-xc)*(gx-xc) + (gy-yc)*(gy-yc);     // offset^2 at bin centre
          // t0 (in samples) of each output sample at offset rng1: 3 fixed-point steps for v(t0)
          for( int j = 1; j <= nt; j++ ) it0[j] = 1 + (int)sqrt( std::max( 0.0, (double)(j-1)*(j-1) - rng1*vel[1] ) );
          for( int iter = 0; iter < 2; iter++ ) {
            for( int j = 1; j <= nt; j++ ) it0[j] = 1 + (int)sqrt( std::max( 0.0, (double)(j-1)*(j-1) - rng1*vel[it0[j]] ) );
          }
          // input sample of each output sample: t_in^2 = t_out^2 + (x_in^2 - x_out^2)/v^2
          for( int j = 1; j <= nt; j++ ) s1[j-1] = (float)( 1.0 + sqrt( std::max( 0.0, (double)(j-1)*(j-1) + (rng2-rng1)*vel[it0[j]] ) ) );
          int n = nt, m = nt;
          spl_( in, &s1[0], &s1[nt], &s1[2*nt], &n, &m );
          memcpy( out, &s1[nt], nt*sizeof(float) );
        }
        else {
          memcpy( out, in, nt*sizeof(float) );
        }
      }
      else {
        ndead++;
        hdrOut->copyFrom( h0 );
        memset( out, 0, nt*sizeof(float) );
      }
      int iiy = ii / ihmt1 + 1;
      int iix = ii % ihmt1 + 1;
      double rang = sqrt( (xc-gxRef)*(xc-gxRef) + (yc-gyRef)*(yc-gyRef) );
      if( vars->idRnum >= 0 ) hdrOut->setDoubleValue( vars->idRnum, iiy );
      if( vars->idTnum >= 0 ) hdrOut->setDoubleValue( vars->idTnum, iix );
      hdrOut->setDoubleValue( vars->idSx, xc );
      hdrOut->setDoubleValue( vars->idSy, yc );
      hdrOut->setDoubleValue( vars->idGx, gxRef );
      hdrOut->setDoubleValue( vars->idGy, gyRef );
      if( vars->idGz >= 0 ) hdrOut->setDoubleValue( vars->idGz, gzRef );
      if( vars->idRang >= 0 )   hdrOut->setDoubleValue( vars->idRang, rang );
      if( vars->idOffset >= 0 ) hdrOut->setDoubleValue( vars->idOffset, rang );
      hdrOut->setDoubleValue( vars->idBinX, xc );
      hdrOut->setDoubleValue( vars->idBinY, yc );
      hdrOut->setIntValue( vars->idRow, iy );
      hdrOut->setIntValue( vars->idCol, ix );
      hdrOut->setFloatValue( vars->idDist, k >= 0 ? (float)rr[ii] : -1.0f );
      hdrOut->setIntValue( vars->idFlag, flag[ii] );
    }
  }

  //--- QC file (replaces COORDS / BCOORDS / BINPLOT)
  if( !vars->qcFile.empty() ) {
    FILE* fq = fopen( vars->qcFile.c_str(), "w" );
    if( fq == NULL ) writer->warning("Cannot open QC file '%s'", vars->qcFile.c_str());
    else {
      fprintf( fq, "# SPK_SHOTBIN QC. flag: 0 dead, 1 trace of the bin, 2 trace from a neighbouring bin\n" );
      fprintf( fq, "#  ix   iy        bin_x        bin_y  hits flag  trace         sou_x         sou_y      dist\n" );
      for( int iy = 1; iy <= ny; iy++ ) {
        for( int ix = 1; ix <= nx; ix++ ) {
          int ii = (iy-1)*nx + ix - 1;
          double xc = xll + (ix-1)*dxx + (iy-1)*dxy;
          double yc = yll + (ix-1)*dyx + (iy-1)*dyy;
          int k = indx[ii];
          fprintf( fq, "%5d %4d %12.2f %12.2f %5d %4d %6d %13.2f %13.2f %9.2f\n", ix, iy, xc, yc, hits[ii], flag[ii],
                   k+1, k >= 0 ? sxv[k] : 0.0, k >= 0 ? syv[k] : 0.0, k >= 0 ? rr[ii] : -1.0 );
        }
      }
      fclose( fq );
      writer->line("  QC file written: %s", vars->qcFile.c_str());
    }
  }

  writer->line("  %d traces available for binning (%d with source outside the grid)", ntrIn, nOutside);
  writer->line("  %d bins created", nbin);
  writer->line("  %d live traces inserted (%d of them from a neighbouring bin)", nlive, nfill);
  writer->line("  %d dead traces inserted", ndead);
  writer->line("  %d duplicate traces (rejected)", ndup);
  if( nlive > 0 ) writer->line("  %.2f average distance from source to bin centre", avgrng/nlive);
  vars->gather->freeAllTraces();
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_spk_shotbin_( csParamDef* pdef ) {
  pdef->setModule( "SPK_SHOTBIN", "Bin the shots of a receiver gather onto a regular grid (SEISPAK XXSHOTBIN)",
                   "Creates a regular rectangular grid from XLL,YLL (lower left corner), the direction to XLR,YLR (lower right corner), "
                   "NX, NY, DX and DY. The trace whose source is closest to a bin centre is moved out to the bin centre with velocity VEL "
                   "(PMO 1) if it is within distance EPS. If no trace is found within EPS, the closest trace of the 8 neighbouring bins is used "
                   "(if within EPS), otherwise a zero trace is placed at the bin. The whole input is ONE receiver gather (OBN/PRM). "
                   "Output: NX*NY traces ordered x fastest, source coordinates = bin centre, headers bin_x, bin_y, row (y), col (x), "
                   "shotbin_dist and shotbin_flag. Parameter names follow SEISPAK pxxshotbin." );

  pdef->addParam( "xll", "X coordinate of lower left corner", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "XLL" );
  pdef->addParam( "yll", "Y coordinate of lower left corner", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "YLL" );
  pdef->addParam( "xlr", "X coordinate of lower right corner (direction of the x axis)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "XLR" );
  pdef->addParam( "ylr", "Y coordinate of lower right corner (direction of the x axis)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "YLR" );
  pdef->addParam( "dx", "Trace spacing in x direction", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "DX" );
  pdef->addParam( "dy", "Trace spacing in y direction", NUM_VALUES_FIXED, "The y axis is 90 degrees counter-clockwise from x" );
  pdef->addValue( "", VALTYPE_NUMBER, "DY" );
  pdef->addParam( "nx", "Number of bins in x direction", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NX" );
  pdef->addParam( "ny", "Number of bins in y direction", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NY" );
  pdef->addParam( "eps", "Tolerance for trace selection to a given bin", NUM_VALUES_FIXED, "Maximum distance source - bin centre" );
  pdef->addValue( "", VALTYPE_NUMBER, "EPS" );

  pdef->addParam( "hmt", "Traces per record of the output", NUM_VALUES_FIXED,
                  "HMT > 1: records of NX traces (one per y-line), else one record of NX*NY traces. "
                  "Sets spk_rnum/spk_tnum if they exist; use the same value in OUTPUT_SEISPAK 'ntraces_record'." );
  pdef->addValue( "1", VALTYPE_NUMBER, "HMT" );

  pdef->addParam( "pmo", "Partial moveout to bin centre", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_OPTION );
  pdef->addOption( "0", "No moveout: the trace is copied to the bin" );
  pdef->addOption( "1", "Hyperbolic moveout from the offset of the trace to the offset of the bin centre, v(t) from VEL" );

  pdef->addParam( "vel", "Velocity-time pair function for moveout to bin centre", NUM_VALUES_VARIABLE,
                  "Pairs v [m/s] and t [ms] (as in SEISPAK: velocity first). Linear interpolation, constant before the first and after the last pair." );
  pdef->addValue( "", VALTYPE_NUMBER, "v1" );
  pdef->addValue( "", VALTYPE_NUMBER, "t1 [ms]" );

  pdef->addParam( "gxbin", "X bin location of the receiver (overrides header value for the moveout)", NUM_VALUES_FIXED, "0: use the receiver header" );
  pdef->addValue( "0", VALTYPE_NUMBER, "GXBIN" );
  pdef->addParam( "gybin", "Y bin location of the receiver (overrides header value for the moveout)", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "GYBIN" );

  pdef->addParam( "hdrs", "Header mnemonics for sorting (source x, y)", NUM_VALUES_FIXED, "Default: spk_sx spk_sy (INPUT_SEISPAK), else sou_x sou_y" );
  pdef->addValue( "", VALTYPE_STRING, "Source x header" );
  pdef->addValue( "", VALTYPE_STRING, "Source y header" );

  pdef->addParam( "hdr_rec", "Receiver x, y, z headers", NUM_VALUES_FIXED, "Default: spk_gx spk_gy spk_gz (or spk_gzel), else rec_x rec_y rec_z" );
  pdef->addValue( "", VALTYPE_STRING, "Receiver x header" );
  pdef->addValue( "", VALTYPE_STRING, "Receiver y header" );
  pdef->addValue( "", VALTYPE_STRING, "Receiver z header" );

  pdef->addParam( "qc_file", "ASCII QC file: one line per bin", NUM_VALUES_FIXED,
                  "ix iy bin_x bin_y, number of sources falling in the bin, flag, input trace, source x y, distance (replaces COORDS/BCOORDS/BINPLOT)" );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
}

//************************************************************************************************
bool start_exec_mod_spk_shotbin_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

void cleanup_mod_spk_shotbin_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  if( vars->gather != NULL ) { delete vars->gather; vars->gather = NULL; }
  delete vars; vars = NULL;
}

extern "C" void _params_mod_spk_shotbin_( csParamDef* pdef ) {
  params_mod_spk_shotbin_( pdef );
}
extern "C" void _init_mod_spk_shotbin_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spk_shotbin_( param, env, writer );
}
extern "C" bool _start_exec_mod_spk_shotbin_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spk_shotbin_( env, writer );
}
extern "C" void _exec_mod_spk_shotbin_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_shotbin_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_shotbin_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spk_shotbin_( env, writer );
}
