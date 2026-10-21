/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/* Module SPK_RESAMPLE: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

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
 * Module: SPK_RESAMPLE
 *
 * SEISPAK program resample.f (and xresample.f): resample every trace to a new
 * sample interval HF with NEWNT samples.
 *  - Spline: cubic spline interpolation, subroutine SPL2 (natural spline whose
 *    tridiagonal system is solved with the fixed coefficients .25, .2667, .2679).
 *  - Decimation: resample.f uses a plain decimation (every 2nd sample) when
 *    HF/SR = 2; here generalised to any integer ratio.
 *  - FFT (method fft): resample2 (G. Hampson) / RESAMP2. The FFT length is the
 *    first even number >= NT with factors 2, 3, 5 only (FAC235).
 *      coarser (HF = R x SR): only the spectrum up to the new Nyquist is kept
 *        and inverse-transformed with length N/R; NEWNT = NT/R as eresample2.f
 *        (eresample2 accepts R = 2**p only; here any integer R).
 *      finer (SR = L x HF): the spectrum is padded with zeros to L x N (RESAMP2);
 *        NEWNT = (NT-1) x L + 1.
 *  - method auto = as resample.f (decimate if HF/SR = 2, otherwise spline);
 *    method spline = as xresample.f.
 * Default NEWNT = (NT-1)*SR/HF + 1, as in the original.
 *
 * SEASEIS changes (marked SEASEIS in the code):
 *  - resample.f passed n = (HF/SR)*NEWNT+1 input samples to SPL2, one or two more
 *    than NT (read past the end of the trace); here n = NT (as xresample.f).
 *  - SPL2 work arrays were fixed (10000/20000 samples); here sized to the trace.
 *  - SPL2 evaluated the polynomial of segment n (not defined) for an output
 *    sample falling exactly on the last input sample; here that sample is y(n).
 *  - RESAMP2 kept only the positive frequencies and took 2 x the real part of the
 *    inverse FFT, which doubles the DC and Nyquist components, and used NT2
 *    without defining it; here the full Hermitian spectrum is used (DC and
 *    Nyquist kept once, Nyquist split in two when padding) and NEWNT is defined.
 *    The execution part of resample2 was not available: the coarser case is the
 *    plain spectral truncation described by eresample2.f.
 *  - Optional anti-alias low-pass (Hann-windowed sinc FIR, zero phase) before a
 *    resampling to a coarser interval. The original has none; default: no.
 */

namespace mod_spk_resample {
  enum { METHOD_AUTO = 0, METHOD_SPLINE = 1, METHOD_DECIMATE = 2, METHOD_FFT = 3 };
  struct VariableStruct {
    int   method;          // resolved: METHOD_SPLINE or METHOD_DECIMATE
    int   ntIn, ntOut;
    float srIn, srOut;     // [ms]
    int   ratio;           // decimation ratio
    // spline set-up (first-call part of SPL2)
    int   m1;
    std::vector<int>    ii;
    std::vector<double> dd;
    double r, r1, r2, r3;
    std::vector<double> y, s, a, b, c;
    // anti-alias
    bool  antialias;
    std::vector<double> fir;    // symmetric, length 2*nh+1
    int   nh;
    std::vector<float>  out;
    // fft
    int nfftIn, nfftOut, nkeep;      // FFT lengths, number of positive bins kept
    bool up;
  };
}
using namespace mod_spk_resample;

//---------------------------------------------------------------- SPL2 (exec part)
static void spl2_exec( VariableStruct* v, float const* yin, float* f ) {
  int n = v->ntIn, m = v->ntOut;
  std::vector<double>& y = v->y; std::vector<double>& s = v->s;
  std::vector<double>& a = v->a; std::vector<double>& b = v->b; std::vector<double>& c = v->c;
  const double b1 = .25, b2 = .2667, b3 = .2679;
  for( int i = 1; i <= n; i++ ) y[i] = yin[i-1];             // 1-based as the original
  // spline coefficients
  s[1] = 0.0;
  for( int i = 2; i <= n-1; i++ ) s[i] = v->r*( y[i+1] - 2.0*y[i] + y[i-1] );
  s[n] = 0.0;
  if( n >= 4 ) {
    s[2] = b1*s[2];
    s[3] = b2*( s[3] - s[2] );
    for( int i = 4; i <= n-1; i++ ) s[i] = b3*( s[i] - s[i-1] );
    for( int i = n-2; i >= 4; i-- ) s[i] = s[i] - b3*s[i+1];
    s[3] = s[3] - b2*s[4];
    s[2] = s[2] - b1*s[3];
  }
  // cubic polynomial coefficients
  for( int j = 1; j <= n-1; j++ ) {
    a[j] = ( s[j+1] - s[j] )*v->r2;
    b[j] = .5*s[j];
    c[j] = ( y[j+1] - y[j] )*v->r1 - ( 2.0*s[j] + s[j+1] )*v->r3;
  }
  // interpolation
  for( int j = 1; j <= v->m1; j++ ) {
    double d1 = v->dd[j];
    int jj = v->ii[j];
    if( jj >= n ) { f[j-1] = (float)y[n]; continue; }       // SEASEIS: last sample
    f[j-1] = (float)( y[jj] + d1*( c[jj] + d1*( b[jj] + d1*a[jj] ) ) );
  }
  for( int j = v->m1+1; j <= m; j++ ) f[j-1] = 0.0f;
}

//---------------------------------------------------------------- FFT (factors 2,3,5)
typedef std::complex<double> cplx;

// FAC235: first even integer >= nx whose only factors are 2, 3 and 5
static int fac235( int nx ) {
  int nxx = ((nx-1)/2)*2;
  for(;;) {
    nxx += 2;
    int n1 = nxx;
    while( n1 % 2 == 0 ) n1 /= 2;
    while( n1 % 3 == 0 ) n1 /= 3;
    while( n1 % 5 == 0 ) n1 /= 5;
    if( n1 == 1 ) return nxx;
  }
}
// unnormalised DFT, sign = -1 forward, +1 inverse; recursive mixed radix (2,3,5)
// w[j] = exp(sign*2*pi*i*j/N) for the full length N; a sub-transform of length n
// uses w[(j*N/n) % N]
static void dft_rec( cplx const* x, int n, int stride, cplx* out, cplx const* w, int N, cplx* tmp ) {
  if( n == 1 ) { out[0] = x[0]; return; }
  int p = ( n % 2 == 0 ) ? 2 : ( n % 3 == 0 ) ? 3 : ( n % 5 == 0 ) ? 5 : n;
  int m = n/p, step = N/n;
  if( p == n ) {                                       // plain DFT (prime length)
    for( int k = 0; k < n; k++ ) {
      cplx acc = 0.0;
      for( int j = 0; j < n; j++ ) acc += x[j*stride]*w[ ((long)j*k*step) % N ];
      out[k] = acc;
    }
    return;
  }
  for( int r = 0; r < p; r++ ) dft_rec( x + r*stride, m, stride*p, tmp + r*m, w, N, out + r*m );
  for( int k = 0; k < m; k++ ) {
    cplx t[5];
    for( int r = 0; r < p; r++ ) t[r] = tmp[r*m + k]*w[ ((long)r*k*step) % N ];   // twiddle
    for( int q = 0; q < p; q++ ) {
      cplx acc = t[0];
      for( int r = 1; r < p; r++ ) acc += t[r]*w[ ((long)r*q*m*step) % N ];        // p-point DFT
      out[k + q*m] = acc;
    }
  }
}
static void dft( std::vector<cplx>& a, int sign ) {
  int N = (int)a.size();
  std::vector<cplx> w( N ), out( N ), tmp( N );
  for( int j = 0; j < N; j++ ) w[j] = std::polar( 1.0, sign*2.0*M_PI*j/(double)N );
  dft_rec( &a[0], N, 1, &out[0], &w[0], N, &tmp[0] );
  a.swap( out );
}

//---------------------------------------------------------------- init
void init_mod_spk_resample_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef* edef = env->execPhaseDef;
  csSuperHeader*  shdr = env->superHeader;
  VariableStruct* v = new VariableStruct();
  edef->setVariables( v );
  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  v->ntIn = shdr->numSamples;
  v->srIn = shdr->sampleInt;

  float hf;
  param->getFloat( "sample_int", &hf );
  if( hf <= 0.0f ) writer->error("sample_int must be > 0");
  v->srOut = hf;

  int newnt = 0;
  if( param->exists("nsamples") ) param->getInt( "nsamples", &newnt );

  string text("auto");
  if( param->exists("method") ) param->getString( "method", &text );
  int method;
  if( text == "auto" ) method = METHOD_AUTO;
  else if( text == "spline" ) method = METHOD_SPLINE;
  else if( text == "decimate" ) method = METHOD_DECIMATE;
  else if( text == "fft" ) method = METHOD_FFT;
  else { writer->error("Unknown option for 'method': %s", text.c_str()); method = METHOD_AUTO; }

  double ratio = hf/v->srIn;
  const double eps = .0001;                                   // as resample.f
  bool isInt = fabs( ratio - floor(ratio + 0.5) ) < eps && ratio >= 1.0 - eps;
  if( method == METHOD_AUTO ) method = ( fabs(ratio - 2.0) < eps ) ? METHOD_DECIMATE : METHOD_SPLINE;
  if( method == METHOD_DECIMATE && !isInt ) {
    writer->error("method decimate needs sample_int = integer x input sample interval (%g / %g = %g)", hf, v->srIn, ratio);
  }
  v->method = method;
  v->ratio  = (int)floor( ratio + 0.5 );

  if( method == METHOD_FFT ) {
    double rr = ( hf < v->srIn ) ? v->srIn/hf : hf/v->srIn;
    if( fabs( rr - floor(rr + 0.5) ) > eps || rr < 1.0 + eps ) {
      writer->error("method fft needs an integer ratio between the sample intervals (%g / %g)", hf, v->srIn);
    }
    int R = (int)floor( rr + 0.5 );
    v->up = ( hf < v->srIn );
    if( v->up ) {                                              // RESAMP2
      v->nfftIn  = fac235( v->ntIn );
      v->nfftOut = R*v->nfftIn;
      if( newnt <= 0 ) newnt = (v->ntIn-1)*R + 1;
    }
    else {                                                     // resample2: N multiple of 2R
      int n = fac235( v->ntIn );
      while( n % (2*R) != 0 ) n = fac235( n+1 );
      v->nfftIn  = n;
      v->nfftOut = n/R;
      if( newnt <= 0 ) newnt = v->ntIn/R;                      // as eresample2.f
    }
    v->nkeep = std::min( v->nfftIn, v->nfftOut )/2;            // bins 0..nkeep (nkeep = Nyquist of the shorter)
    if( newnt > v->nfftOut ) writer->error("nsamples %d larger than the FFT output length %d", newnt, v->nfftOut);
  }
  if( newnt <= 0 ) newnt = (int)( (v->ntIn-1)*v->srIn/hf + 1 );    // as resample.f
  v->ntOut = newnt;

  // anti-alias
  text = "no";
  if( param->exists("antialias") ) param->getString( "antialias", &text );
  v->antialias = ( text == "yes" );
  if( text != "yes" && text != "no" ) writer->error("Unknown option for 'antialias': %s", text.c_str());
  v->nh = 0;
  if( v->antialias ) {
    if( ratio <= 1.0 + eps ) {
      writer->warning("antialias ignored: output sample interval is not coarser than the input");
      v->antialias = false;
    }
    else {
      double fnyqOut = 500.0/hf;                               // [Hz]
      double fc = 0.8*fnyqOut;
      if( param->exists("aa_freq") ) { float x; param->getFloat( "aa_freq", &x ); fc = x; }
      double fnyqIn = 500.0/v->srIn;
      if( fc <= 0.0 || fc >= fnyqIn ) writer->error("aa_freq must be between 0 and the input Nyquist (%g Hz)", fnyqIn);
      double wc = fc/fnyqIn;                                   // cut-off / input Nyquist
      v->nh = (int)ceil( 8.0*ratio );
      v->fir.assign( 2*v->nh+1, 0.0 );
      double sum = 0.0;
      for( int k = -v->nh; k <= v->nh; k++ ) {
        double x = M_PI*wc*k;
        double sinc = ( k == 0 ) ? 1.0 : sin(x)/x;
        double win  = 0.5*( 1.0 + cos( M_PI*k/(double)(v->nh+1) ) );
        v->fir[k+v->nh] = wc*sinc*win; sum += v->fir[k+v->nh];
      }
      for( size_t k = 0; k < v->fir.size(); k++ ) v->fir[k] /= sum;  // unit gain at 0 Hz
      writer->line("  Anti-alias: low-pass %.1f Hz (output Nyquist %.1f Hz), FIR %d taps", fc, fnyqOut, 2*v->nh+1);
    }
  }

  // spline set-up (first-call part of SPL2)
  if( v->method == METHOD_SPLINE ) {
    int n = v->ntIn, m = v->ntOut;                             // SEASEIS: n = NT
    double r = hf/v->srIn, hy = v->srIn;
    v->m1 = std::min( m, 1 + (int)( (n-1)/r ) );
    v->ii.assign( v->m1+1, 0 ); v->dd.assign( v->m1+1, 0.0 );
    for( int j = 1; j <= v->m1; j++ ) {
      double x = 1.0 + r*(j-1);
      v->ii[j] = (int)x;
      v->dd[j] = hy*( x - v->ii[j] );
    }
    v->r  = 6.0/(hy*hy);
    v->r1 = 1.0/hy;
    v->r2 = 1.0/(6.0*hy);
    v->r3 = hy/6.0;
    v->y.assign( n+2, 0.0 ); v->s.assign( n+2, 0.0 );
    v->a.assign( n+2, 0.0 ); v->b.assign( n+2, 0.0 ); v->c.assign( n+2, 0.0 );
  }
  v->out.assign( std::max( v->ntOut, v->ntIn ), 0.0f );

  writer->line("  Input : %d samples, %g ms (%g ms)", v->ntIn, v->srIn, (v->ntIn-1)*v->srIn);
  writer->line("  Output: %d samples, %g ms (%g ms)", v->ntOut, v->srOut, (v->ntOut-1)*v->srOut);
  writer->line("  Method: %s", v->method == METHOD_SPLINE ? "cubic spline (SPL2)" :
               v->method == METHOD_FFT ? "FFT (resample2 / RESAMP2)" : "decimation");
  if( v->method == METHOD_FFT ) writer->line("  FFT length in/out %d / %d (FAC235)", v->nfftIn, v->nfftOut);
  if( v->method == METHOD_DECIMATE ) writer->line("  Decimation ratio %d (keeps 1 of every %d samples)", v->ratio, v->ratio);
  if( v->method == METHOD_FFT && v->antialias ) {
    writer->line("  antialias ignored: method fft already removes the spectrum above the new Nyquist");
    v->antialias = false;
  }
  if( ratio > 1.0 + eps && !v->antialias && v->method != METHOD_FFT )
    writer->line("  No anti-alias filter (as SEISPAK): energy above %.1f Hz is aliased", 500.0/hf);

  shdr->numSamples = v->ntOut;
  shdr->sampleInt  = v->srOut;
}

//---------------------------------------------------------------- exec
void exec_mod_spk_resample_( csTraceGather* traceGather, int* port, int* numTrcToKeep,
                             csExecPhaseEnv* env, csLogWriter* writer )
{
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  float* smp = traceGather->trace(0)->getTraceSamples();
  int n = v->ntIn;

  std::vector<float> filt;
  float const* in = smp;
  if( v->antialias ) {                                         // SEASEIS: optional low-pass
    filt.assign( n, 0.0f );
    for( int i = 0; i < n; i++ ) {
      double acc = 0.0;
      for( int k = -v->nh; k <= v->nh; k++ ) {
        int j = i + k;
        if( j < 0 || j >= n ) continue;
        acc += v->fir[k+v->nh]*smp[j];
      }
      filt[i] = (float)acc;
    }
    in = &filt[0];
  }

  float* out = &v->out[0];
  if( v->method == METHOD_DECIMATE ) {
    for( int k = 0; k < v->ntOut; k++ ) {
      int j = k*v->ratio;                                      // tr2(k) = tr1(ratio*(k-1)+1)
      out[k] = ( j < n ) ? in[j] : 0.0f;                       // SEASEIS: zero past the end
    }
  }
  else if( v->method == METHOD_FFT ) {
    int N = v->nfftIn, M = v->nfftOut, K = v->nkeep;
    std::vector<cplx> X( N, 0.0 );
    for( int i = 0; i < n && i < N; i++ ) X[i] = in[i];
    dft( X, -1 );
    std::vector<cplx> Y( M, 0.0 );
    Y[0] = X[0];
    for( int k = 1; k < K; k++ ) { Y[k] = X[k]; Y[M-k] = X[N-k]; }
    if( v->up ) { Y[K] = 0.5*X[K]; Y[M-K] = 0.5*X[K]; }       // input Nyquist split in two
    else        { Y[K] = cplx( X[K].real(), 0.0 ); }           // new Nyquist (real)
    dft( Y, +1 );
    double fac = 1.0/(double)N;                                // keeps amplitudes
    for( int k = 0; k < v->ntOut; k++ ) out[k] = (float)( fac*Y[k].real() );
  }
  else {
    spl2_exec( v, in, out );
  }
  memcpy( smp, out, v->ntOut*sizeof(float) );
}

//---------------------------------------------------------------- params
void params_mod_spk_resample_( csParamDef* pdef ) {
  pdef->setModule( "SPK_RESAMPLE", "Resample traces (SEISPAK resample / xresample)",
                   "Cubic spline interpolation (SEISPAK SPL2) or plain decimation to a new sample interval. "
                   "method auto = as resample.f (decimation if the new interval is twice the input one, else spline). "
                   "SEISPAK applies no anti-alias filter; use antialias yes for a low-pass before a coarser resampling." );

  pdef->addParam( "sample_int", "Output sample interval [ms] (HF)", NUM_VALUES_FIXED );
  pdef->addValue( "", VALTYPE_NUMBER, "Sample interval [ms]" );

  pdef->addParam( "nsamples", "Number of output samples (NEWNT)", NUM_VALUES_FIXED,
                  "Default: (NT-1)*SR/HF + 1 (same trace length); method fft: NT/R (coarser) or (NT-1)*L+1 (finer)" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Number of samples (0 = default)" );

  pdef->addParam( "method", "Resampling method", NUM_VALUES_FIXED );
  pdef->addValue( "auto", VALTYPE_OPTION );
  pdef->addOption( "auto", "As resample.f: decimation if HF = 2 x SR, otherwise cubic spline" );
  pdef->addOption( "spline", "Cubic spline (SPL2) always, as xresample.f" );
  pdef->addOption( "decimate", "Keep every N-th sample (HF = N x SR, N integer)" );
  pdef->addOption( "fft", "FFT resampling (resample2 / RESAMP2): integer ratio, coarser or finer" );

  pdef->addParam( "antialias", "Anti-alias low-pass before a coarser resampling (not in SEISPAK)", NUM_VALUES_FIXED );
  pdef->addValue( "no", VALTYPE_OPTION );
  pdef->addOption( "no", "No filter (as SEISPAK)" );
  pdef->addOption( "yes", "Zero-phase windowed-sinc low-pass" );

  pdef->addParam( "aa_freq", "Anti-alias cut-off frequency [Hz]", NUM_VALUES_FIXED,
                  "Default: 0.8 x output Nyquist" );
  pdef->addValue( "", VALTYPE_NUMBER, "Frequency [Hz]" );
}

bool start_exec_mod_spk_resample_( csExecPhaseEnv* env, csLogWriter* writer ) { return true; }

void cleanup_mod_spk_resample_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* v = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  delete v; v = NULL;
}

extern "C" void _params_mod_spk_resample_( csParamDef* pdef ) { params_mod_spk_resample_( pdef ); }
extern "C" void _init_mod_spk_resample_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_spk_resample_( param, env, writer );
}
extern "C" bool _start_exec_mod_spk_resample_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_spk_resample_( env, writer );
}
extern "C" void _exec_mod_spk_resample_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_spk_resample_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_spk_resample_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_spk_resample_( env, writer );
}
