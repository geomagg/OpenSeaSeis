/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module SPK_MUTE: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

#include "cseis_includes.h"
#include <cstring>
#include <cmath>
#include <vector>
#include <complex>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

/**
 * CSEIS - Seabed Seismic Processing System
 * Module: SPK_MUTE
 *
 * SEISPAK process MUTE (emute.f / mute.f, with XTINTERP, JOEHINT, SPLYN, XSHIFT1, CLP):
 *  1. Hyperbolic mute: t = sqrt( (T0+(k-1)*DT0)**2 + (RANG*1000/VEL)**2 ) [ms]
 *     RANG = offset header (IX0 = 0) or sqrt( ((i-IX0)*DX)**2 + ((k-IY0)*DY)**2 )
 *     (k = record, i = trace in record). T0 >= 0: top mute (zero above t, 10-sample
 *     linear ramp up to t). T0 < 0: bottom mute (zero from t to the end).
 *  2. SHIFT [ms]: whole-sample shift, or FFT phase shift for a fractional one;
 *     the samples shifted in through the edge are zeroed. LRM: output end time.
 *  3. XTM mutes: per record block, each pair of (trace, time) curves defines a
 *     region (below the first curve and above the second) that is zeroed.
 *     Curve interpolation: 0 linear, -1 cubic spline, -2 smoothed linear.
 *  4. CLIP (subroutine CLP): CLIP < 0 zeroes a +-100 ms window around every
 *     sample with |amp| > |CLIP|; CLIP > 0 clips the amplitudes to +-CLIP.
 *
 * SEASEIS changes (marked SEASEIS in the code):
 *  - CLIP was read but CLP was never called; it is applied here (last step).
 *  - A shift smaller than one sample was skipped (test on the integer part only).
 *  - The fractional shift used only the positive frequencies and 2 x the real part
 *    of the inverse FFT, doubling the DC and Nyquist components; here the full
 *    Hermitian spectrum is shifted.
 *  - The whole-sample shift filled the vacated samples with the first sample of
 *    the trace; here with zeros (as the fractional shift).
 *  - When LRM made the output longer than the input, the extra samples kept old
 *    buffer contents; here they are zero.
 *  - XTM sample indices are limited to the trace length.
 *  - Headers GY, GZ/GZEL and MUTS were read but not used: not used here either.
 * Record and trace numbers come from trace headers (default spk_rec, spk_trc as
 * written by INPUT_SEISPAK). For single-trace records (HMT = 1) SEISPAK uses the
 * record number as trace number: set hdr_trc accordingly.
 */

namespace mod_spk_mute {
  struct XtmPair {
    int ixl, ixr;                       // traces (1-based) where both curves exist
    std::vector<double> tt1, tt2;       // times [ms] indexed by trace (1-based)
  };
  struct XtmBlock {
    int rec1, rec2;
    std::vector<XtmPair> pairs;
  };
  struct VariableStruct {
    int    lrm, ix0, iy0;
    double shift, t0, dx, dy, vel, dt0, clip;
    int    it0sign;                     // +1 top mute, -1 bottom mute
    int    nt, ntt;
    float  sr;
    int    idRec, idTrc, idOff;
    int    ishft; double ashft;
    bool   doShift, fracShift;
    int    nfft;
    std::vector< std::complex<double> > phase;   // exp(-i w ashft), full spectrum
    std::vector<XtmBlock> blocks;
    int    lwin2;
    std::vector<float> work;
  };
}
using namespace mod_spk_mute;
typedef std::complex<double> cplx;

//---------------------------------------------------------------- FFT (power of 2, radix 2)
static void fft2( std::vector<cplx>& a, int sign ) {
  int n = (int)a.size();
  for( int i = 1, j = 0; i < n; i++ ) {
    int bit = n >> 1;
    for( ; j & bit; bit >>= 1 ) j ^= bit;
    j ^= bit;
    if( i < j ) std::swap( a[i], a[j] );
  }
  for( int len = 2; len <= n; len <<= 1 ) {
    cplx wl = std::polar( 1.0, sign*2.0*M_PI/len );
    for( int i = 0; i < n; i += len ) {
      cplx w = 1.0;
      for( int k = 0; k < len/2; k++ ) {
        cplx u = a[i+k], v = a[i+k+len/2]*w;
        a[i+k] = u + v; a[i+k+len/2] = u - v; w *= wl;
      }
    }
  }
}

//---------------------------------------------------------------- JOEHINT (1-based arrays)
// linear (mode 0, -3) or smoothed linear (mode -2) interpolation
static void joehint( double const* x, double const* y, double const* xx, double* yy, int mm, int nn, int mode ) {
  std::vector<double> dy( mm+2, 0.0 );
  int mm1 = mm - 1;
  int mmm = mm1;
  if( mode == -2 ) mmm = 1;
  for( int i = 1; i <= mmm; i++ ) dy[i] = ( y[i+1] - y[i] )/( x[i+1] - x[i] );
  dy[mm] = ( y[mm] - y[mm-1] )/( x[mm] - x[mm-1] );
  if( !( mm < 3 || mode == -3 || mode == 0 ) ) {
    for( int i = 2; i <= mm1; i++ ) {                  // slopes at interior points
      double x0 = x[i], xm = x[i-1]-x0, xp = x[i+1]-x0;
      double fm = y[i-1], f0 = y[i], fp = y[i+1];
      dy[i] = 0.0;
      double pm = fm - f0, pp = fp - f0;
      double wm = xm*xm + pm*pm;
      double xm2 = xm*xm/wm, xm3 = xm*xm2, xm4 = xm*xm3, xm5 = xm*xm4;
      double xp2 = xp*xp/wm, xp3 = xp*xp2, xp4 = xp*xp3, xp5 = xp*xp4;   // (wm, as the original)
      double a1 = .25*(xm4-xp4), b1 = .333333*(xm3-xp3), g1 = .333333*(pm*xm2-pp*xp2);
      double a2 = .20*(xm5-xp5), b2 = a1, g2 = .25*(pm*xm3-pp*xp3);
      double den = a2*b1 - a1*b2, top = a2*g1 - a1*g2;
      if( !( den*1000000. <= top ) ) dy[i] = top/den;
    }
  }
  int j = 1, ii = 1;
  double a1 = 0, a2 = 0, a3 = 0, a4 = 0, eps = 0;
  for( int i = 1; i <= nn; i++ ) {
    for(;;) {
      if( ii == i ) {
        double delx = x[j+1] - x[j];
        eps = .000001*delx;
        double dxr = 1.0/(delx*delx);
        a1 = y[j]; a2 = dy[j];
        a3 = ( (y[j+1]-a1) - (x[j+1]-x[j])*a2 )*dxr;
        a4 = ( (dy[j+1]-a2) - 2.*(x[j+1]-x[j])*a3 )*dxr;
        if( mode != -2 ) { a3 = 0.; a4 = 0.; }
      }
      double del1 = xx[i] - x[j], del2 = x[j+1] - xx[i];
      if( ( del1+eps >= 0. && del2+eps >= 0. ) || ( del1 < 0. && j == 1 ) || j == mm1 ) {
        yy[i] = a1 + del1*( a2 + del1*( a3 - del2*a4 ) );
        break;
      }
      j = j + 1; ii = i;
    }
  }
}

//---------------------------------------------------------------- SPLYN (1-based arrays)
static void splyn( double const* x, double const* xx, double const* y, double* f, int n, int m ) {
  int n1 = n-1, n2 = n-2, n3 = n-3;
  if( n1 <= 1 ) {
    double den = 1.0/( x[2] - x[1] );
    for( int i = 1; i <= m; i++ ) {
      double cl = ( x[2] - xx[i] )*den, cr = ( xx[i] - x[1] )*den;
      f[i] = cl*y[1] + cr*y[2];
    }
    return;
  }
  std::vector<double> a( n+2, 0.0 ), b( n+2, 0.0 ), c( n+2, 0.0 ), d( n+2, 0.0 );
  for( int i = 1; i <= n2; i++ ) {
    double h1 = x[i+1]-x[i], h2 = x[i+2]-x[i+1];
    a[i] = h1*h1*h2; c[i] = h1*h2*h2; b[i] = 2.*( a[i] + c[i] );
    d[i+1] = 6.*( h1*( y[i+2]-y[i+1] ) - h2*( y[i+1]-y[i] ) );
  }
  double c6 = 1./6.;
  b[1] = 1./b[1]; c[1] = c[1]*b[1]; d[2] = d[2]*b[1];
  if( n3 != 0 ) {
    if( n3 != 1 ) {
      for( int i = 2; i <= n3; i++ ) {
        b[i] = 1./( b[i] - c[i-1]*a[i-1] );
        c[i] = b[i]*c[i];
        d[i+1] = b[i]*( d[i+1] - a[i-1]*d[i] );
      }
    }
    b[n2] = b[n2] - c[n3]*a[n3];
    d[n2+1] = ( d[n2+1] - a[n3]*d[n2] )/b[n2];
    for( int k = 1; k <= n3; k++ ) { int i = n2-k; d[i+1] = d[i+1] - c[i]*d[i+2]; }
  }
  d[1] = 0.; d[n] = 0.;
  for( int i = 1; i <= n1; i++ ) {
    double hi = x[i+1]-x[i], hib = 1./hi;
    b[i] = .5*d[i]; a[i] = ( d[i+1]-d[i] )*c6*hib;
    c[i] = hib*( y[i+1]-y[i] ) - c6*( hi*( 2.*d[i] + d[i+1] ) );
  }
  double eps = ( x[2]-x[1] )*.0001;
  int j = 1;
  for( int i = 1; i <= m; i++ ) {
    for(;;) {
      double del1 = xx[i]-x[j], del2 = x[j+1]-xx[i];
      if( ( del1+eps >= 0. && del2+eps >= 0. ) || ( del1 < 0. && j == 1 ) || j == n1 ) {
        f[i] = y[j] + del1*( c[j] + del1*( b[j] + del1*a[j] ) );
        break;
      }
      j++;
    }
  }
}

//---------------------------------------------------------------- XTM parsing (as mute.f)
static double xv( std::vector<double> const& xtm, int j ) {          // 1-based, zero past the end
  return ( j >= 1 && j <= (int)xtm.size() ) ? xtm[j-1] : 0.0;
}
// reads one event descriptor starting at jxtm (x t x t ... mode); fills tt over its traces
static bool xt_event( std::vector<double> const& xtm, int& jxtm, std::vector<double>& tt, int& ixl, int& ixr, csLogWriter* writer ) {
  std::vector<double> x( 1, 0.0 ), t( 1, 0.0 );
  do {
    x.push_back( xv(xtm, jxtm) ); t.push_back( xv(xtm, jxtm+1) );
    jxtm += 2;
    if( jxtm > (int)xtm.size() + 2 ) return false;
  } while( xv(xtm, jxtm) > 0.0 );
  int mode = (int)xv( xtm, jxtm );
  jxtm += 1;
  int np = (int)x.size() - 1;
  ixl = (int)( x[1] + .1 );
  ixr = (int)( x[np] + .0001 );
  int nnx = ixr - ixl + 1;
  if( np < 2 || nnx < 1 ) { writer->warning("XTM event with less than 2 points or decreasing trace numbers (trace %g): ignored", x[1]); return true; }
  std::vector<double> xx( nnx+1 ), yy( nnx+1 );
  for( int i = 1; i <= nnx; i++ ) xx[i] = ixl + i - 1;
  if( mode <= -2 || mode == 0 ) joehint( &x[0], &t[0], &xx[0], &yy[0], np, nnx, mode );
  else splyn( &x[0], &xx[0], &t[0], &yy[0], np, nnx );
  if( (int)tt.size() < ixr+1 ) tt.resize( ixr+1, tt.empty() ? 0.0 : tt.back() );
  for( int i = 1; i <= nnx; i++ ) if( ixl+i-1 >= 1 ) tt[ixl+i-1] = yy[i];
  writer->line("    curve: %d points, traces %d-%d, %s", np, ixl, ixr, mode == 0 ? "linear" : mode == -2 ? "smoothed linear" : "spline");
  return true;
}

//---------------------------------------------------------------- init
void init_mod_spk_mute_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  csTraceHeaderDef* hdef = env->headerDef;
  VariableStruct* v = new VariableStruct();
  edef->setVariables( v );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  // defaults as emute.f
  v->lrm = 0; v->clip = 0.; v->ix0 = 0; v->iy0 = 1; v->dx = 0.; v->dy = 0.;
  v->t0 = 0.; v->dt0 = 0.; v->vel = 1000000.; v->shift = 0.;
  float f;
  if( param->exists("lrm") )   { param->getFloat("lrm", &f);   v->lrm = (int)f; }
  if( param->exists("t0") )    { param->getFloat("t0", &f);    v->t0 = f; }
  if( param->exists("dt0") )   { param->getFloat("dt0", &f);   v->dt0 = f; }
  if( param->exists("ix0") )   { param->getFloat("ix0", &f);   v->ix0 = (int)f; }
  if( param->exists("iy0") )   { param->getFloat("iy0", &f);   v->iy0 = (int)f; }
  if( param->exists("dx") )    { param->getFloat("dx", &f);    v->dx = f; }
  if( param->exists("dy") )    { param->getFloat("dy", &f);    v->dy = f; }
  if( param->exists("shift") ) { param->getFloat("shift", &f); v->shift = f; }
  if( param->exists("clip") )  { param->getFloat("clip", &f);  v->clip = f; }
  if( param->exists("vel") )   { param->getFloat("vel", &f);   v->vel = f; }
  if( v->vel == 0.0 ) writer->error("VEL must not be zero");

  v->nt = shdr->numSamples;
  v->sr = shdr->sampleInt;
  double sr = v->sr;
  v->ntt = v->nt;
  if( v->lrm > 0 ) v->ntt = (int)( (double)v->lrm/sr + 1.1 );
  if( v->t0 < 0. ) { v->t0 = -v->t0; v->it0sign = -1; } else v->it0sign = 1;

  // headers
  string name;
  name = "spk_rec"; if( param->exists("hdr_rec") ) param->getString("hdr_rec", &name);
  if( !hdef->headerExists(name) ) writer->error("Trace header '%s' (record number) not found: set hdr_rec", name.c_str());
  v->idRec = hdef->headerIndex( name );
  name = "spk_trc"; if( param->exists("hdr_trc") ) param->getString("hdr_trc", &name);
  if( !hdef->headerExists(name) ) writer->error("Trace header '%s' (trace number) not found: set hdr_trc", name.c_str());
  v->idTrc = hdef->headerIndex( name );
  v->idOff = -1;
  if( v->ix0 == 0 ) {
    if( param->exists("hdr_offset") ) param->getString("hdr_offset", &name);
    else name = hdef->headerExists("spk_rang") ? "spk_rang" : "offset";
    if( !hdef->headerExists(name) ) writer->error("Offset header '%s' not found (IX0 = 0 uses the offset header): set hdr_offset", name.c_str());
    v->idOff = hdef->headerIndex( name );
  }

  // shift
  v->ashft = v->shift/sr;
  v->ishft = (int)v->ashft;                                  // truncation, as Fortran
  double rshft = v->ashft - v->ishft;
  if( fabs(rshft) < .001 ) rshft = 0.;
  v->doShift   = ( v->ishft != 0 || rshft != 0. );            // SEASEIS: shift < 1 sample applied
  v->fracShift = ( rshft != 0. );
  if( v->fracShift ) {
    int log2 = (int)( log( (double)v->nt - 1. )/log( 2. ) ) + 1;
    v->nfft = 1 << log2;
    v->phase.assign( v->nfft, 1.0 );
    for( int k = 0; k < v->nfft; k++ ) {
      int kk = ( k <= v->nfft/2 ) ? k : k - v->nfft;          // signed frequency index
      double w = 2.0*M_PI*kk/v->nfft;
      v->phase[k] = std::polar( 1.0, -w*v->ashft );            // delay by ashft samples
      if( k == v->nfft/2 ) v->phase[k] = cos( w*v->ashft );    // Nyquist: real part
    }
  }

  // XTM
  std::vector<double> xtm;
  int nlines = param->getNumLines("xtm");
  for( int il = 0; il < nlines; il++ ) {
    csVector<string> vals; param->getAll( "xtm", &vals, il );
    for( int k = 0; k < vals.size(); k++ ) xtm.push_back( atof( vals.at(k).c_str() ) );
  }
  if( !xtm.empty() ) {
    writer->line("  XTM mutes (%d values):", (int)xtm.size());
    int jxtm = 1;
    while( xv(xtm, jxtm) > 0.0 ) {                            // record block
      XtmBlock blk;
      blk.rec1 = (int)xv(xtm, jxtm); jxtm++;
      if( xv(xtm, jxtm) < 0.0 ) { blk.rec2 = (int)(-xv(xtm, jxtm)); jxtm++; } else blk.rec2 = blk.rec1;
      writer->line("   records %d-%d", blk.rec1, blk.rec2);
      do {                                                    // pairs of events
        XtmPair pr;
        std::vector<double> tt1, tt2;
        int l1, r1, l2, r2;
        if( !xt_event( xtm, jxtm, tt1, l1, r1, writer ) || !xt_event( xtm, jxtm, tt2, l2, r2, writer ) ) {
          writer->error("XTM list ends in the middle of an event (see the format in the module help)");
        }
        pr.ixl = std::max( 1, std::max( l1, l2 ) );
        pr.ixr = std::min( r1, r2 );
        pr.tt1 = tt1; pr.tt2 = tt2;
        if( pr.ixl <= pr.ixr ) blk.pairs.push_back( pr );
        writer->line("    -> zero between the two curves on traces %d-%d", pr.ixl, pr.ixr);
      } while( xv(xtm, jxtm) > 0.0 );
      jxtm++;                                                 // block terminator 0
      v->blocks.push_back( blk );
    }
  }

  v->lwin2 = (int)( 100./sr + 1 );
  v->work.assign( std::max( v->nt, v->ntt ) + 8, 0.0f );

  writer->line("  T0 %g ms (%s)  DT0 %g  VEL %g  IX0 %d IY0 %d DX %g DY %g", v->t0, v->it0sign > 0 ? "top mute" : "bottom mute",
               v->dt0, v->vel, v->ix0, v->iy0, v->dx, v->dy);
  writer->line("  SHIFT %g ms (%.3f samples%s)  LRM %d ms -> %d output samples", v->shift, v->ashft,
               v->fracShift ? ", FFT phase shift" : "", v->lrm, v->ntt);
  if( v->clip < 0 ) writer->line("  CLIP %g: zero +-%d samples around |amp| > %g", v->clip, v->lwin2, -v->clip);
  else if( v->clip > 0 ) writer->line("  CLIP %g: amplitudes clipped to +-%g", v->clip, v->clip);
  shdr->numSamples = v->ntt;
}

//---------------------------------------------------------------- exec
void exec_mod_spk_mute_( csTraceGather* traceGather, int* port, int* numTrcToKeep,
                         csExecPhaseEnv* env, csLogWriter* writer )
{
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csTrace* trace = traceGather->trace(0);
  csTraceHeader* th = trace->getTraceHeader();
  float* s = trace->getTraceSamples();
  int nt = v->nt, ntt = v->ntt;
  double sr = v->sr;
  int k = th->intValue( v->idRec );
  int i = th->intValue( v->idTrc );

  //--- 1. hyperbolic mute
  double rang;
  if( v->ix0 == 0 ) rang = th->doubleValue( v->idOff );
  else rang = sqrt( pow( (i - v->ix0)*v->dx, 2 ) + pow( (k - v->iy0)*v->dy, 2 ) );
  double x = pow( v->t0 + (k-1)*v->dt0, 2 ) + pow( rang*1000./v->vel, 2 );
  int it = (int)( 1.0 + sqrt(x)/sr );
  it = std::min( it, nt );
  if( it > 1 ) {
    if( v->it0sign == 1 ) {
      for( int j = 1; j <= it; j++ ) {
        float scale = std::max( 0.0f, 1.0f + (float)(j - it)/10.0f );
        s[j-1] *= scale;
      }
    }
    else {
      for( int j = it; j <= nt; j++ ) s[j-1] = 0.0f;
    }
  }

  //--- 2. shift / output length
  float* w = &v->work[0];
  for( int j = 0; j < (int)v->work.size(); j++ ) w[j] = 0.0f;
  if( !v->doShift ) {
    for( int j = 0; j < std::min( nt, ntt ); j++ ) w[j] = s[j];       // SEASEIS: zero past NT
  }
  else if( !v->fracShift ) {
    int sh = v->ishft;
    for( int j = 0; j < nt; j++ ) {                                    // SEASEIS: zero fill
      int jd = j + sh;
      if( jd >= 0 && jd < ntt ) w[jd] = s[j];
    }
  }
  else {
    int n = v->nfft;
    std::vector<cplx> y( n, 0.0 );
    for( int j = 0; j < nt; j++ ) y[j] = s[j];
    fft2( y, -1 );
    for( int kk = 0; kk < n; kk++ ) y[kk] *= v->phase[kk];
    fft2( y, +1 );
    for( int j = 0; j < std::min( nt, ntt ); j++ ) w[j] = (float)( y[j].real()/n );
    if( v->ishft > 0 ) { for( int j = 1; j <= v->ishft; j++ ) w[j-1] = 0.0f; }
    else { for( int j = nt + v->ishft; j <= nt; j++ ) if( j >= 1 && j <= ntt ) w[j-1] = 0.0f; }
  }
  memcpy( s, w, ntt*sizeof(float) );

  //--- 3. XTM mutes
  for( size_t b = 0; b < v->blocks.size(); b++ ) {
    XtmBlock const& blk = v->blocks[b];
    if( k < blk.rec1 || k > blk.rec2 ) continue;
    for( size_t p = 0; p < blk.pairs.size(); p++ ) {
      XtmPair const& pr = blk.pairs[p];
      if( i < pr.ixl || i > pr.ixr ) continue;
      int it1 = (int)( pr.tt1[i]/sr + 1.5 );
      int it2 = (int)( pr.tt2[i]/sr + 1.5 );
      if( it1 < it2 ) {
        it1 = std::max( it1, 1 ); it2 = std::min( it2, ntt );         // SEASEIS: limits
        for( int j = it1; j <= it2; j++ ) s[j-1] = 0.0f;
      }
    }
  }

  //--- 4. CLIP (CLP)                                                 // SEASEIS: CLP applied
  if( v->clip != 0.0 ) {
    double aclip = fabs( v->clip );
    for( int ii = 1; ii <= ntt; ii++ ) {
      if( fabs( s[ii-1] ) > aclip && v->clip < 0. ) {
        int j1 = std::max( 1, ii - v->lwin2 ), j2 = std::min( ntt, ii + v->lwin2 );
        for( int j = j1; j <= j2; j++ ) s[j-1] = 0.0f;
      }
      if( v->clip > 0. ) {
        if( s[ii-1] >  aclip ) s[ii-1] = (float)aclip;
        if( s[ii-1] < -aclip ) s[ii-1] = (float)(-aclip);
      }
    }
  }
}

//---------------------------------------------------------------- params
void params_mod_spk_mute_( csParamDef* pdef ) {
  pdef->setModule( "SPK_MUTE", "Mute, shift and clip (SEISPAK mute)",
                   "Hyperbolic mute t = sqrt((T0+(k-1)*DT0)^2 + (RANG*1000/VEL)^2) [ms] (T0 > 0 top mute with 10-sample ramp, "
                   "T0 < 0 bottom mute), static SHIFT, output end time LRM, XTM mutes between pairs of trace-time curves, and CLIP. "
                   "k = record number (hdr_rec), traces numbered by hdr_trc." );
  pdef->addParam( "t0", "T0 for the mute [ms] (< 0: bottom mute)", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "T0 [ms]" );
  pdef->addParam( "dt0", "Record-to-record change in T0 [ms]", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "DT0 [ms]" );
  pdef->addParam( "vel", "Moveout velocity for the mute [m/s]", NUM_VALUES_FIXED );
  pdef->addValue( "1000000", VALTYPE_NUMBER, "VEL [m/s]" );
  pdef->addParam( "ix0", "Trace number of the T0 location (0: use the offset header)", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "IX0" );
  pdef->addParam( "iy0", "Line (record) number of the T0 location", NUM_VALUES_FIXED );
  pdef->addValue( "1", VALTYPE_NUMBER, "IY0" );
  pdef->addParam( "dx", "Trace-to-trace distance for the moveout [m]", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "DX [m]" );
  pdef->addParam( "dy", "Record-to-record (line) distance for the moveout [m]", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "DY [m]" );
  pdef->addParam( "shift", "Shift applied to every trace [ms] (> 0: later)", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "SHIFT [ms]" );
  pdef->addParam( "lrm", "End time of the output section [ms] (0: unchanged)", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "LRM [ms]" );
  pdef->addParam( "clip", "CLIP < 0: zero +-100 ms around |amp| > |CLIP|; CLIP > 0: clip to +-CLIP", NUM_VALUES_FIXED );
  pdef->addValue( "0", VALTYPE_NUMBER, "CLIP (0 = off)" );
  pdef->addParam( "xtm", "X-T mutes, SEISPAK format; may be given on several lines", NUM_VALUES_VARIABLE,
                  "Blocks: record [-last_record], then pairs of curves 'trace time trace time ... code' "
                  "(code 0 linear, -1 spline, -2 smoothed linear), block ended by 0. "
                  "The region below the first curve and above the second is zeroed. Example: "
                  "xtm 1 -5 10 500 20 500 30 600 -2 10 600 30 800 0 0" );
  pdef->addValue( "", VALTYPE_NUMBER, "Values" );
  pdef->addParam( "hdr_rec", "Trace header with the record number (k)", NUM_VALUES_FIXED );
  pdef->addValue( "spk_rec", VALTYPE_STRING, "Header name" );
  pdef->addParam( "hdr_trc", "Trace header with the trace number in the record (i)", NUM_VALUES_FIXED );
  pdef->addValue( "spk_trc", VALTYPE_STRING, "Header name" );
  pdef->addParam( "hdr_offset", "Trace header with the offset (SEISPAK RANG), used when IX0 = 0", NUM_VALUES_FIXED,
                  "Default: spk_rang if present, otherwise offset" );
  pdef->addValue( "spk_rang", VALTYPE_STRING, "Header name" );
}

bool start_exec_mod_spk_mute_( csExecPhaseEnv* env, csLogWriter* writer ) { return true; }

void cleanup_mod_spk_mute_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  delete v; v = NULL;
}

extern "C" void _params_mod_spk_mute_( csParamDef* pdef ) { params_mod_spk_mute_( pdef ); }
extern "C" void _init_mod_spk_mute_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spk_mute_( param, env, writer );
}
extern "C" bool _start_exec_mod_spk_mute_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spk_mute_( env, writer );
}
extern "C" void _exec_mod_spk_mute_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_mute_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_mute_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spk_mute_( env, writer );
}
