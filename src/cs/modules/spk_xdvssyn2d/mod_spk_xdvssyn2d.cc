/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module SPK_XDVSSYN2D: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

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
 * Module: SPK_XDVSSYN2D
 *
 * Input module: 2D acoustic time-domain modelling of one shot from SEISPAK program
 * XDVSSYN2DRS (syn3d): pseudo-spectral / finite-difference Laplacian, optional variable
 * density, KUTTA terms of the cos(dt*sqrt(-L)) expansion per time step and Koslov absorbing
 * boundaries. Output: section recorded at depth IDP, wavefield snapshots or the velocity model.
 */

extern "C" {
  void spkrs_( float* vel, float* rho, int* iden, int* nx, int* nz, int* nt, float* dt, float* dx, float* dz,
               float* apx, int* kutta, float* raa, int* nx0, int* nz0, int* idp, float* gam1, float* fs,
               int* itype, int* itzr, float* wave, int* ntw, float* out, int* nout, int* ierr );
}

namespace mod_spk_xdvssyn2d {
  struct VariableStruct {
    int nx, nz, nt, nx0, nz0, idp, kutta, itzr, itype, iden;
    float sr, dx, ddz, apx, raa, gam1, fs;
    int nsnap;
    int nOutTraces, nOutSamples;
    std::vector<float> vel, rho, wave, out;
    int traceCounter;
    bool atEOF, computed;
    int hdrId_trcno, hdrId_ffid, hdrId_chan, hdrId_soux, hdrId_souz, hdrId_recx, hdrId_recz, hdrId_offset, hdrId_time;
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
    else return "only NBYTE 2 and 4 supported for model files";
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

  // Model on the grid (nz,nx) as velcom4/rhocom4: trace ix of the first record, sample iz
  // (limited to ntt-ilim). SEASEIS: a 2D file with one trace per record uses the records as x.
  void loadModel( std::string const& fname, char const* what, int ilim, int nx, int nz, float ddz,
                  std::vector<float>* model, csLogWriter* writer ) {
    int hmt, hmr, ntt; float smp;
    std::vector<float> dat;
    std::string err = readSeispak( fname, &hmt, &hmr, &ntt, &smp, &dat );
    if( !err.empty() ) {
      std::string msg = std::string(what) + " file '" + fname + "': " + err;
      writer->error( "%s", msg.c_str() );
    }
    int ntrx = hmt;
    bool recAsX = ( hmt == 1 && hmr > 1 );
    if( recAsX ) ntrx = hmr;
    writer->line("  %s file: %s  (%d traces x %d records, %d samples, sample interval %g)%s", what, fname.c_str(), hmt, hmr, ntt, smp,
                 recAsX ? "  -> one trace per record: records used as x" : "");
    if( ntrx < nx ) writer->warning("%s file has %d locations in x, model has NX = %d: last location repeated", what, ntrx, nx);
    if( ntt - ilim < nz ) writer->warning("%s file has %d samples, model has NZ = %d: last sample repeated", what, ntt, nz);
    if( fabs( smp - ddz ) > 0.01f*ddz ) writer->warning("%s file sample interval (%g) differs from DDZ (%g). Samples are used as grid points (as in SEISPAK): resample the file to DDZ.", what, smp, ddz);
    model->assign( (size_t)nx*nz, 0.0f );
    int nlim = std::max( 1, ntt - ilim );
    for( int ix = 0; ix < nx; ix++ ) {
      int it = std::min( ix+1, ntrx ) - 1;
      for( int iz = 0; iz < nz; iz++ ) {
        int izz = std::min( iz+1, nlim ) - 1;
        (*model)[(size_t)ix*nz + iz] = dat[(size_t)it*ntt + izz];
      }
    }
  }
}
using namespace mod_spk_xdvssyn2d;

//*************************************************************************************************
// Init phase
//*************************************************************************************************
void init_mod_spk_xdvssyn2d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
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
  vars->computed = false;
  vars->nx = vars->nz = vars->nt = 0;
  vars->nx0 = vars->nz0 = vars->idp = 0;
  vars->sr = 0.0f; vars->dx = 0.0f; vars->ddz = 0.0f;
  vars->apx = 0.0f; vars->kutta = 2; vars->raa = 0.0f; vars->gam1 = 20.0f; vars->fs = 0.0f;
  vars->itzr = 10; vars->itype = 2; vars->iden = 0; vars->nsnap = 0;

  std::string text;
  if( param->exists("nx") )    param->getInt( "nx", &vars->nx );
  if( param->exists("nz") )    param->getInt( "nz", &vars->nz );
  if( param->exists("nt") )    param->getInt( "nt", &vars->nt );
  if( param->exists("sr") )    param->getFloat( "sr", &vars->sr );
  if( param->exists("dx") )    param->getFloat( "dx", &vars->dx );
  if( param->exists("ddz") )   param->getFloat( "ddz", &vars->ddz );
  if( param->exists("nx0") )   param->getInt( "nx0", &vars->nx0 );
  if( param->exists("nz0") )   param->getInt( "nz0", &vars->nz0 );
  if( param->exists("idp") )   param->getInt( "idp", &vars->idp );
  if( param->exists("apx") )   param->getFloat( "apx", &vars->apx );
  if( param->exists("kutta") ) param->getInt( "kutta", &vars->kutta );
  if( param->exists("raa") )   param->getFloat( "raa", &vars->raa );
  if( param->exists("gam1") )  param->getFloat( "gam1", &vars->gam1 );
  if( param->exists("itzr") )  param->getInt( "itzr", &vars->itzr );
  if( param->exists("free_surface") ) {
    param->getString( "free_surface", &text );
    if( !text.compare("yes") ) vars->fs = 1.0f;
    else if( text.compare("no") ) writer->error("Unknown option for 'free_surface': %s", text.c_str());
  }
  if( param->exists("output") ) {
    param->getString( "output", &text );
    if( !text.compare("section") ) vars->itype = 2;
    else if( !text.compare("snapshots") ) vars->itype = 1;
    else if( !text.compare("velocity") ) vars->itype = 5;
    else writer->error("Unknown option for 'output': %s", text.c_str());
  }

  if( vars->nx <= 0 )  writer->error("Required parameter NX (grid points in x) not given");
  if( vars->nz <= 0 )  writer->error("Required parameter NZ (grid points in z) not given");
  if( vars->nt <= 0 )  writer->error("Required parameter NT (time steps) not given");
  if( vars->sr <= 0 )  writer->error("Required parameter SR (time step [ms]) not given");
  if( vars->dx <= 0 )  writer->error("Required parameter DX not given");
  if( vars->ddz <= 0 ) writer->error("Required parameter DDZ not given");
  if( vars->nx % 2 != 0 || vars->nz % 2 != 0 ) writer->error("NX and NZ must be even (the Fortran operators work on pairs of traces)");
  if( vars->nx < 8 || vars->nz < 8 ) writer->error("NX and NZ must be >= 8");
  if( vars->nx0 < 1 || vars->nx0 > vars->nx ) writer->error("Source position NX0 (1..NX) required");
  if( vars->nz0 < 1 || vars->nz0 > vars->nz ) writer->error("Source depth NZ0 (1..NZ) required");
  if( vars->itype == 2 && (vars->idp < 1 || vars->idp > vars->nz) ) writer->error("Receiver depth IDP (1..NZ) required");
  float a = vars->apx;
  if( !( a == 0.0f || a == 1.0f || a == 2.0f || a == 3.0f || (a >= 4.0f && a <= 7.0f && a == floorf(a)) ) )
    writer->error("APX = %g not supported. Use 0/2 (pseudo-spectral), 1 (2D FFT), 3 (staggered with density), 4..7 (FD in z, half-length APX)", a);
  if( vars->kutta < 1 ) writer->error("KUTTA must be >= 1 (terms of the time expansion)");
  int nn = ( (int)lrintf(vars->gam1) == 0 ) ? 20 : (int)lrintf(fabsf(vars->gam1));
  if( nn > 50 ) writer->error("GAM1 (absorbing boundary width) must be <= 50");
  if( 2*nn >= vars->nx || 2*(int)fabsf(vars->gam1)+nn >= vars->nz ) writer->error("Grid too small for the absorbing boundaries (width %d)", nn);
  if( vars->raa == 1.0f && ( vars->nz0 < 3 || vars->nz0 > vars->nz-2 || vars->nx0 < 2 || vars->nx0 > vars->nx-1 ) )
    writer->error("RAA 1 needs 3 <= NZ0 <= NZ-2 and 2 <= NX0 <= NX-1");
  if( vars->itype == 1 && vars->itzr < 1 ) writer->error("ITZR (snapshot interval) must be >= 1");

  //---------------------------------------------
  // Velocity, density, wavelet
  if( !param->exists("vel_file") ) writer->error("Velocity required: 'vel_file' (SEISPAK file with v(z) traces, sampled at DDZ, as ONLVELS)");
  std::string fname;
  param->getString( "vel_file", &fname );
  loadModel( fname, "Velocity", 2, vars->nx, vars->nz, vars->ddz, &vars->vel, writer );
  if( param->exists("rho_file") ) {
    param->getString( "rho_file", &fname );
    loadModel( fname, "Density", 0, vars->nx, vars->nz, vars->ddz, &vars->rho, writer );
    vars->iden = 1;
    for( size_t i = 0; i < vars->rho.size(); i++ ) if( !(vars->rho[i] > 0.0f) ) writer->error("Density must be > 0");
  }
  else {
    vars->rho.assign( 1, 1.0f );
  }
  float vmin = 1e30f, vmax = 0.0f;
  for( size_t i = 0; i < vars->vel.size(); i++ ) {
    if( !(vars->vel[i] > 0.0f) ) writer->error("Velocity must be > 0");
    vmin = std::min( vmin, vars->vel[i] ); vmax = std::max( vmax, vars->vel[i] );
  }
  if( param->exists("wavelet_file") ) {
    param->getString( "wavelet_file", &fname );
    FILE* fp = fopen( fname.c_str(), "r" );
    if( fp == NULL ) writer->error("Cannot open wavelet file '%s'", fname.c_str());
    float w;
    while( fscanf( fp, "%f", &w ) == 1 ) vars->wave.push_back( w );
    fclose( fp );
    if( vars->wave.size() < 2 ) writer->error("Wavelet file '%s': need at least 2 samples", fname.c_str());
    writer->line("  Wavelet file: %s  (%d samples at SR)", fname.c_str(), (int)vars->wave.size());
  }
  if( vars->wave.empty() ) vars->wave.assign( 1, 0.0f );

  //---------------------------------------------
  // Stability of the time stepping: 2*cos(x) approximated by KUTTA terms of its series,
  // x = v*dt*|k|max (pseudo-spectral), must stay within [-2,2]
  float dt = 0.001f * vars->sr;
  double xmax = vmax * dt * M_PI * sqrt( 1.0/(vars->dx*vars->dx) + 1.0/(vars->ddz*vars->ddz) );
  bool stable = true;
  for( int i = 1; i <= 200 && stable; i++ ) {
    double x = xmax * i / 200.0, term = 2.0, sum = 2.0;
    for( int k = 1; k <= vars->kutta; k++ ) { term *= -x*x / ( (2.0*k)*(2.0*k-1.0) ); sum += term; }
    if( sum < -2.0 - 1e-9 || sum > 2.0 + 1e-9 ) stable = false;
  }
  double ppw = vmin / ( 60.0 * std::max( vars->dx, vars->ddz ) );   // points per wavelength at 60 Hz
  writer->line("  Velocity %g - %g   v*dt*kmax = %.3f (KUTTA %d)%s", vmin, vmax, xmax, vars->kutta, stable ? "" : "  UNSTABLE");
  if( !stable ) writer->warning("Time stepping is unstable (v*dt*kmax = %.3f too large for KUTTA %d): reduce SR or increase KUTTA", xmax, vars->kutta);
  if( ppw < 2.5 ) writer->warning("Only %.1f grid points per wavelength at 60 Hz for the lowest velocity: expect dispersion", ppw);

  //---------------------------------------------
  shdr->domain = DOMAIN_XT;
  if( vars->itype == 2 ) {
    vars->nOutTraces = vars->nx; vars->nOutSamples = vars->nt;
    shdr->numSamples = vars->nt;
    shdr->sampleInt  = vars->sr;
  }
  else if( vars->itype == 1 ) {
    vars->nsnap = (vars->nt-1)/vars->itzr + 1;
    vars->nOutTraces = vars->nx * vars->nsnap; vars->nOutSamples = vars->nz;
    shdr->numSamples = vars->nz;
    shdr->sampleInt  = vars->ddz;
    shdr->domain = DOMAIN_XD;
  }
  else {
    vars->nOutTraces = vars->nx; vars->nOutSamples = vars->nz;
    shdr->numSamples = vars->nz;
    shdr->sampleInt  = vars->ddz;
    shdr->domain = DOMAIN_XD;
  }
  vars->hdrId_trcno  = hdef->addStandardHeader( HDR_TRCNO.name );
  vars->hdrId_ffid   = hdef->addStandardHeader( HDR_FFID.name );
  vars->hdrId_chan   = hdef->addStandardHeader( HDR_CHAN.name );
  vars->hdrId_soux   = hdef->addStandardHeader( HDR_SOU_X.name );
  vars->hdrId_souz   = hdef->addStandardHeader( HDR_SOU_Z.name );
  vars->hdrId_recx   = hdef->addStandardHeader( HDR_REC_X.name );
  vars->hdrId_recz   = hdef->addStandardHeader( HDR_REC_Z.name );
  vars->hdrId_offset = hdef->addStandardHeader( HDR_OFFSET.name );
  vars->hdrId_time = -1;
  if( vars->itype == 1 ) {
    if( !hdef->headerExists("spk_time") ) hdef->addHeader( TYPE_FLOAT, "spk_time", "Time of the wavefield snapshot [ms]" );
    vars->hdrId_time = hdef->headerIndex( "spk_time" );
  }

  char const* apxText = ( a == 1.0f ) ? "2D FFT" : ( a == 3.0f && vars->iden ) ? "staggered pseudo-spectral" :
                        ( a >= 4.0f && !vars->iden ) ? "x pseudo-spectral, z finite difference" : "pseudo-spectral";
  writer->line("  Grid NX %d  NZ %d  DX %g  DDZ %g   NT %d  SR %g ms   APX %g (%s)  KUTTA %d  %s",
               vars->nx, vars->nz, vars->dx, vars->ddz, vars->nt, vars->sr, a, apxText, vars->kutta,
               vars->iden ? "variable density" : "constant density");
  writer->line("  Source NX0 %d NZ0 %d (x %g, z %g)  RAA %g   Receivers at IDP %d (z %g)   GAM1 %g  free surface %s",
               vars->nx0, vars->nz0, (vars->nx0-1)*vars->dx, (vars->nz0-1)*vars->ddz, vars->raa, vars->idp,
               (vars->idp-1)*vars->ddz, vars->gam1, vars->fs != 0.0f ? "yes" : "no");
  if( vars->itype == 1 ) writer->line("  Output: %d snapshots (every %d steps), each NX traces x NZ samples", vars->nsnap, vars->itzr);
  else if( vars->itype == 5 ) writer->line("  Output: velocity model used (no modelling)");
  else writer->line("  Output: shot record, NX traces x NT samples (first sample = time SR)");
}

//*************************************************************************************************
// Exec phase
//*************************************************************************************************
void exec_mod_spk_xdvssyn2d_(
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
  int ntr = vars->nOutTraces, ns = vars->nOutSamples;
  if( !vars->computed ) {
    vars->out.assign( (size_t)ntr * ns, 0.0f );
    int nout = (int)vars->out.size();
    int ierr = 0;
    float dt = 0.001f * vars->sr;
    int ntw = (int)vars->wave.size();
    if( ntw < 2 ) ntw = 0;
    writer->line("  Computing %d time steps...", vars->nt);
    spkrs_( &vars->vel[0], &vars->rho[0], &vars->iden, &vars->nx, &vars->nz, &vars->nt, &dt, &vars->dx, &vars->ddz,
            &vars->apx, &vars->kutta, &vars->raa, &vars->nx0, &vars->nz0, &vars->idp, &vars->gam1, &vars->fs,
            &vars->itype, &vars->itzr, &vars->wave[0], &ntw, &vars->out[0], &nout, &ierr );
    if( ierr == 2 ) writer->error("Not enough memory for the modelling grid");
    float amax = 0.0f;
    for( size_t i = 0; i < vars->out.size(); i++ ) {
      if( !std::isfinite( vars->out[i] ) ) writer->error("Modelling diverged (NaN/Inf): reduce SR or increase KUTTA");
      amax = std::max( amax, fabsf( vars->out[i] ) );
    }
    writer->line("  Done. Max amplitude %g", amax);
    vars->vel.clear(); vars->rho.clear();
    vars->computed = true;
  }
  int k = vars->traceCounter;
  if( k == ntr-1 ) vars->atEOF = true;
  csTrace* trace = traceGather->trace(0);
  memcpy( trace->getTraceSamples(), &vars->out[(size_t)k*ns], ns*sizeof(float) );
  csTraceHeader* trcHdr = trace->getTraceHeader();
  int ix = k % vars->nx + 1;
  int irec = k / vars->nx + 1;
  double soux = (vars->nx0-1)*vars->dx, recx = (ix-1)*vars->dx;
  trcHdr->setIntValue( vars->hdrId_trcno, k+1 );
  trcHdr->setIntValue( vars->hdrId_ffid, irec );
  trcHdr->setIntValue( vars->hdrId_chan, ix );
  trcHdr->setDoubleValue( vars->hdrId_soux, soux );
  trcHdr->setFloatValue( vars->hdrId_souz, (vars->nz0-1)*vars->ddz );
  trcHdr->setDoubleValue( vars->hdrId_recx, recx );
  trcHdr->setFloatValue( vars->hdrId_recz, vars->idp > 0 ? (vars->idp-1)*vars->ddz : 0.0f );
  trcHdr->setFloatValue( vars->hdrId_offset, (float)(recx - soux) );
  if( vars->hdrId_time >= 0 ) trcHdr->setFloatValue( vars->hdrId_time, (1 + (irec-1)*vars->itzr) * vars->sr );
  vars->traceCounter += 1;
}

//********************************************************************************
// Parameter definition
//********************************************************************************
void params_mod_spk_xdvssyn2d_( csParamDef* pdef ) {
  pdef->setModule( "SPK_XDVSSYN2D", "2D acoustic time-domain shot modelling (SEISPAK XDVSSYN2DRS)",
                   "Input module. Propagates one shot in a 2D velocity (and optional density) grid with the XDVSSYN2DRS "
                   "operators: pseudo-spectral or finite-difference Laplacian, KUTTA terms of the time expansion and Koslov "
                   "absorbing boundaries. The source is an initial wavefield at (NX0,NZ0) (RAA) or an injected wavelet. "
                   "Output: shot record at depth IDP, wavefield snapshots, or the velocity grid. "
                   "Computation is done by the original Fortran routines." );
  pdef->addParam( "nx", "Grid points in x (even)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NX" );
  pdef->addParam( "nz", "Grid points in z (even)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NZ" );
  pdef->addParam( "dx", "Grid spacing in x", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "DX" );
  pdef->addParam( "ddz", "Grid spacing in z", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "DDZ" );
  pdef->addParam( "nt", "Number of time steps (= output samples)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NT" );
  pdef->addParam( "sr", "Time step [ms] (DT of SEISPAK in seconds x 1000)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "SR" );
  pdef->addParam( "nx0", "Source position in x (grid point, 1 = first)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NX0" );
  pdef->addParam( "nz0", "Source position in z (grid point, 1 = top)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "NZ0" );
  pdef->addParam( "idp", "Receiver depth (grid point) for the shot record", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "IDP" );
  pdef->addParam( "apx", "Spatial operator", NUM_VALUES_FIXED, "With rho_file: 0,2,4..7 test2dr, 1 rhsfftr, 3 staggered grid (xtest3drs)" );
  pdef->addValue( "0", VALTYPE_NUMBER, "0/2/3: pseudo-spectral (xrhsfft), 1: 2D FFT Laplacian (rhsfft), 4..7: x pseudo-spectral and z finite difference of half-length APX (Holberg)" );
  pdef->addParam( "kutta", "Terms of the time expansion of cos(dt*sqrt(-L))", NUM_VALUES_FIXED, "1 = 2nd order in time; more terms allow a larger time step" );
  pdef->addValue( "2", VALTYPE_NUMBER, "KUTTA" );
  pdef->addParam( "raa", "Source type (initial wavefield)", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "0: smoothed spike, 1: vertical dipole pattern, 2: spike with k-space taper, 3/4: derivative of Gaussian (SPGAUSSr/u), "
                  "negative: unfiltered spike (and -3/-4 Gaussian)" );
  pdef->addParam( "gam1", "Width of the absorbing boundaries (grid points)", NUM_VALUES_FIXED, "0 = 20. The bottom 2*GAM1 points get the surface velocity" );
  pdef->addValue( "20", VALTYPE_NUMBER, "GAM1" );
  pdef->addParam( "free_surface", "Top boundary", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "Absorbing boundary at the top too (as SEISPAK with FS = 0)" );
  pdef->addOption( "yes", "No absorbing taper at the top (FS <> 0)" );
  pdef->addParam( "output", "What to output", NUM_VALUES_FIXED );
  pdef->addValue( "section", VALTYPE_OPTION );
  pdef->addOption( "section", "Shot record at depth IDP: NX traces, NT samples (ITYPE 2)" );
  pdef->addOption( "snapshots", "Wavefield every ITZR steps: one ensemble (ffid) of NX traces x NZ depth samples per snapshot (ITYPE 1)" );
  pdef->addOption( "velocity", "Velocity grid used (ITYPE 5)" );
  pdef->addParam( "itzr", "Time steps between snapshots", NUM_VALUES_FIXED );
  pdef->addValue( "10", VALTYPE_NUMBER, "ITZR" );
  pdef->addParam( "vel_file", "SEISPAK file with velocity v(z) traces (ONLVELS)", NUM_VALUES_FIXED,
                  "Trace ix of the first record = grid column ix, sample iz = grid point iz (the file must be sampled at DDZ). "
                  "A 2D file with one trace per record (e.g. from SPK_VXT2ONL + OUTPUT_SEISPAK) uses the records as x." );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addParam( "rho_file", "SEISPAK file with density (ONLRHOS), same layout as vel_file", NUM_VALUES_FIXED, "If not given: constant density" );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
  pdef->addParam( "wavelet_file", "Text file with wavelet samples at SR, injected at the source (file WAVELET of SEISPAK)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_STRING, "File name" );
}

//************************************************************************************************
bool start_exec_mod_spk_xdvssyn2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}
void cleanup_mod_spk_xdvssyn2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  delete vars; vars = NULL;
}
extern "C" void _params_mod_spk_xdvssyn2d_( csParamDef* pdef ) {
  params_mod_spk_xdvssyn2d_( pdef );
}
extern "C" void _init_mod_spk_xdvssyn2d_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spk_xdvssyn2d_( param, env, writer );
}
extern "C" bool _start_exec_mod_spk_xdvssyn2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spk_xdvssyn2d_( env, writer );
}
extern "C" void _exec_mod_spk_xdvssyn2d_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_xdvssyn2d_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_xdvssyn2d_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spk_xdvssyn2d_( env, writer );
}
