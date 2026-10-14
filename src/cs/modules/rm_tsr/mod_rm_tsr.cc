/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */
/*
 * Module RM_TSR: Two-Sided Recursive (TSR) inverse filter of a known wavelet.
 * Author: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026.
 * Module prefix RM = Roland Marschall.
 *
 * Convention: z = unit delay, B(z) = sum_k b_k z^k (wavelet b_0..b_n, sampled at the trace interval).
 *
 *   TSR(z) = 1/B(z) = B(1/z) / [ B(z) B(1/z) ] = B(1/z) / R(z) = B(1/z) / [ B0(z) B0(1/z) ]
 *
 *   R(z)  : autocorrelation of b (+ white noise on r_0)
 *   B0(z) : minimum-phase factor of R, obtained by Wiener-Levinson (spike filter R a = e0, a = 1/(b00 B0))
 *
 * The TSR filter is the unique STABLE inverse of B, whether or not B is minimum phase:
 *   1) v_t = sum_k b_k x_{t+k}                        correlation with b (anticausal FIR, carries the phase)
 *   2) u_t = ( v_t - sum_{k>=1} b0_k u_{t-k} ) / b0_0  causal recursion,      t = 0 -> N
 *   3) y_t = ( u_t - sum_{k>=1} b0_k y_{t+k} ) / b0_0  anticausal recursion,  t = N -> 0
 * Steps 2-3 can alternatively be applied with the truncated FIR 1/B0 from Levinson (method fir).
 *
 * Modes:  tsr       full inverse 1/B (amplitude 1/|B|, phase -phi_B)
 *         whiten    two-sided zero-phase filter 1/R = 1/|B|^2 only (no correlation step)
 *         minphase  causal 1/B0 only (classic spike deconvolution with the known wavelet; wrong phase if B is not min. phase)
 */

#include "cseis_includes.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using namespace cseis_system;
using namespace cseis_geolib;
using namespace std;

namespace mod_rm_tsr {
  static const int MODE_TSR      = 1;
  static const int MODE_WHITEN   = 2;
  static const int MODE_MINPHASE = 3;
  static const int METHOD_RECURSIVE = 11;
  static const int METHOD_FIR       = 12;

  struct VariableStruct {
    int mode;
    int method;
    int nb;                 // wavelet length (n+1)
    int nInv;               // length of FIR 1/B0 (Wiener/Levinson order)
    int shift;              // output shift [samples] = zero_time/dt (mode tsr only)
    int nPad;               // zero padding on each side [samples]
    double* b;              // wavelet b_k, k=0..nb-1
    double* b0;             // minimum-phase factor B0, k=0..nb-1
    double* inv;            // FIR 1/B0, k=0..nInv-1
    double* buf1;           // work buffers (padded trace)
    double* buf2;
    int nBuf;
    int numTraces;
  };

  //--------------------------------------------------------------------
  // Levinson recursion: solve the Toeplitz system R f = (E,0,...,0)^T with f_0 = 1 (prediction-error filter).
  // r[k], k=0..nr-1 (r_k = 0 for k >= nr). Returns E. f has length m.
  double levinson( double const* r, int nr, int m, double* f ) {
    std::vector<double> tmp( m, 0.0 );
    for( int i = 0; i < m; i++ ) f[i] = 0.0;
    f[0] = 1.0;
    double E = r[0];
    for( int k = 1; k < m; k++ ) {
      double acc = 0.0;
      for( int j = 0; j < k; j++ ) {
        int lag = k - j;
        if( lag < nr ) acc += f[j] * r[lag];
      }
      double c = -acc / E;
      for( int j = 0; j <= k; j++ ) tmp[j] = f[j] + c * f[k-j];
      for( int j = 0; j <= k; j++ ) f[j] = tmp[j];
      E *= ( 1.0 - c*c );
      if( E <= 0.0 ) return E;
    }
    return E;
  }
}
using namespace mod_rm_tsr;

//*************************************************************************************************
// Init phase
//
void init_mod_rm_tsr_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer )
{
  csExecPhaseDef*   edef = env->execPhaseDef;
  csSuperHeader*    shdr = env->superHeader;
  VariableStruct* vars   = new VariableStruct();
  edef->setVariables( vars );

  edef->setTraceSelectionMode( TRCMODE_FIXED, 1 );

  vars->mode   = MODE_TSR;
  vars->method = METHOD_RECURSIVE;
  vars->nb = vars->nInv = vars->shift = vars->nPad = vars->nBuf = 0;
  vars->b = vars->b0 = vars->inv = vars->buf1 = vars->buf2 = NULL;
  vars->numTraces = 0;

  if( shdr->domain != DOMAIN_XT ) {
    writer->error("RM_TSR requires time-domain (XT) input.");
  }
  double dt = shdr->sampleInt;

  //--- Wavelet: parameter list and/or ASCII file
  std::vector<double> wav;
  if( param->exists("wavelet") && param->exists("filename") ) {
    writer->error("Specify either 'wavelet' or 'filename', not both.");
  }
  if( param->exists("wavelet") ) {
    int nLines = param->getNumLines("wavelet");
    for( int iline = 0; iline < nLines; iline++ ) {
      int nv = param->getNumValues( "wavelet", iline );
      for( int iv = 0; iv < nv; iv++ ) {
        double val;
        param->getDoubleAtLine( "wavelet", &val, iline, iv );
        wav.push_back( val );
      }
    }
  }
  else if( param->exists("filename") ) {
    std::string filename;
    param->getString( "filename", &filename );
    FILE* fin = fopen( filename.c_str(), "r" );
    if( fin == NULL ) writer->error("Cannot open wavelet file '%s'.", filename.c_str());
    char line[4096];
    while( fgets( line, 4096, fin ) != NULL ) {
      char* hash = strchr( line, '#' );
      if( hash != NULL ) *hash = '\0';
      char* tok = strtok( line, " \t,;\r\n" );
      while( tok != NULL ) {
        char* end;
        double val = strtod( tok, &end );
        if( end == tok ) {
          fclose( fin );
          writer->error("Wavelet file '%s': cannot read number '%s'.", filename.c_str(), tok);
        }
        wav.push_back( val );
        tok = strtok( NULL, " \t,;\r\n" );
      }
    }
    fclose( fin );
  }
  else {
    writer->error("Wavelet not specified. Use parameter 'wavelet' (list of samples) or 'filename' (ASCII file).");
  }
  if( wav.size() < 2 ) writer->error("Wavelet must have at least 2 samples (found %d).", (int)wav.size());
  if( wav[0] == 0.0 ) {
    writer->warning("First wavelet sample is zero: leading zeros only delay the wavelet. Consider removing them and using 'zero_time'.");
  }

  //--- Options
  if( param->exists("mode") ) {
    std::string text;
    param->getString( "mode", &text );
    text = toLowerCase( text );
    if( !text.compare("tsr") ) vars->mode = MODE_TSR;
    else if( !text.compare("whiten") ) vars->mode = MODE_WHITEN;
    else if( !text.compare("minphase") ) vars->mode = MODE_MINPHASE;
    else writer->error("Unknown option for parameter 'mode': '%s'", text.c_str());
  }
  if( param->exists("method") ) {
    std::string text;
    param->getString( "method", &text );
    text = toLowerCase( text );
    if( !text.compare("recursive") ) vars->method = METHOD_RECURSIVE;
    else if( !text.compare("fir") ) vars->method = METHOD_FIR;
    else writer->error("Unknown option for parameter 'method': '%s'", text.c_str());
  }
  double whiteNoise = 0.1;
  if( param->exists("white_noise") ) param->getDouble( "white_noise", &whiteNoise );
  if( whiteNoise < 0.0 ) writer->error("White noise must be >= 0: %f %%", whiteNoise);

  int nb = (int)wav.size();
  int nInv = 4 * nb;
  if( param->exists("filter_length") ) {
    double len;
    param->getDouble( "filter_length", &len );
    nInv = (int)floor( len / dt + 0.5 ) + 1;
    if( nInv < nb ) writer->error("Filter length (%f ms = %d samples) must be at least the wavelet length (%d samples).", len, nInv, nb);
  }
  double zeroTime = 0.0;
  if( param->exists("zero_time") ) param->getDouble( "zero_time", &zeroTime );
  int shift = (int)floor( zeroTime / dt + 0.5 );
  if( shift < 0 || shift >= nb ) writer->error("Zero time (%f ms = sample %d) must be inside the wavelet (0 ... %d samples).", zeroTime, shift, nb-1);
  if( vars->mode != MODE_TSR ) shift = 0;

  vars->nb    = nb;
  vars->nInv  = nInv;
  vars->shift = shift;
  vars->b     = new double[nb];
  vars->b0    = new double[nb];
  vars->inv   = new double[nInv];
  for( int k = 0; k < nb; k++ ) vars->b[k] = wav[k];

  //--- Autocorrelation R (+ white noise)
  std::vector<double> r( nb, 0.0 );
  for( int k = 0; k < nb; k++ ) {
    double sum = 0.0;
    for( int i = 0; i+k < nb; i++ ) sum += wav[i] * wav[i+k];
    r[k] = sum;
  }
  if( r[0] <= 0.0 ) writer->error("Wavelet has zero energy.");
  double r0orig = r[0];
  r[0] *= ( 1.0 + whiteNoise/100.0 );

  //--- Wiener-Levinson: prediction-error filter f (f_0=1), R f = (E,0..0)  =>  1/B0 = f/sqrt(E),  B0_0 = sqrt(E)
  std::vector<double> f( nInv );
  double E = levinson( &r[0], nb, nInv, &f[0] );
  if( E <= 0.0 ) writer->error("Levinson recursion failed (E = %g). Increase 'white_noise'.", E);
  double sqE = sqrt( E );
  for( int k = 0; k < nInv; k++ ) vars->inv[k] = f[k] / sqE;

  //--- B0 = inverse of (1/B0), truncated to nb coefficients
  vars->b0[0] = sqE;
  for( int t = 1; t < nb; t++ ) {
    double sum = 0.0;
    for( int k = 1; k <= t && k < nInv; k++ ) sum += f[k] * vars->b0[t-k];
    vars->b0[t] = -sum;
  }

  //--- Quality checks
  // (1) factorization error: autocorrelation of B0 vs R
  double errR = 0.0, normR = 0.0;
  for( int k = 0; k < nb; k++ ) {
    double sum = 0.0;
    for( int i = 0; i+k < nb; i++ ) sum += vars->b0[i] * vars->b0[i+k];
    errR  += (sum - r[k]) * (sum - r[k]);
    normR += r[k] * r[k];
  }
  errR = sqrt( errR / normR );
  // (2) tail energy of FIR 1/B0 (last 10%): truncation indicator
  double eTot = 0.0, eTail = 0.0;
  int iTail = nInv - nInv/10 - 1;
  for( int k = 0; k < nInv; k++ ) {
    eTot += vars->inv[k]*vars->inv[k];
    if( k >= iTail ) eTail += vars->inv[k]*vars->inv[k];
  }
  double tailRatio = eTail / eTot;
  // (3) recursion 1/B0 reproduces FIR 1/B0, and decays (stability)
  int nTest = 4 * nInv;
  std::vector<double> h( nTest, 0.0 );
  for( int t = 0; t < nTest; t++ ) {
    double sum = (t == 0) ? 1.0 : 0.0;
    for( int k = 1; k < nb && k <= t; k++ ) sum -= vars->b0[k] * h[t-k];
    h[t] = sum / vars->b0[0];
  }
  double errH = 0.0;
  for( int k = 0; k < nInv; k++ ) errH += (h[k]-vars->inv[k])*(h[k]-vars->inv[k]);
  errH = sqrt( errH / eTot );
  double eLast = 0.0;
  for( int t = 3*nInv; t < nTest; t++ ) eLast += h[t]*h[t];
  // Recursion is accepted if B0 factorizes R and its impulse response decays (B0 minimum phase).
  // errH (recursion vs truncated FIR) is reported only: it is dominated by the FIR truncation.
  bool recursionOk = ( errR < 1.0e-3 ) && std::isfinite( eLast ) && ( eLast < 1.0e-6 * eTot );
  // (4) is b itself minimum phase? (b == B0 up to white noise)
  double diffMin = 0.0, normB = 0.0;
  for( int k = 0; k < nb; k++ ) {
    diffMin += (wav[k]-vars->b0[k])*(wav[k]-vars->b0[k]);
    normB   += wav[k]*wav[k];
  }
  diffMin = sqrt( diffMin / normB );

  if( vars->method == METHOD_RECURSIVE && !recursionOk ) {
    writer->warning("Recursive 1/B0 is not reliable (factorization error %.2e, or impulse response does not decay). Switching to method 'fir'. Consider a longer 'filter_length' or more 'white_noise'.", errR);
    vars->method = METHOD_FIR;
  }
  if( tailRatio > 1.0e-4 ) {
    writer->warning("FIR 1/B0 has not decayed at the end of the filter (tail energy %.2e). Wavelet roots are close to the unit circle: increase 'filter_length' and/or 'white_noise'.", tailRatio);
  }
  if( errR > 1.0e-3 ) {
    writer->warning("Spectral factorization error |B0 B0* - R| / |R| = %.2e. Increase 'filter_length'.", errR);
  }

  //--- Padding: covers the decay of the anticausal/causal parts and the correlation length
  vars->nPad = nInv + nb + shift;
  vars->nBuf = shdr->numSamples + 2 * vars->nPad;
  vars->buf1 = new double[vars->nBuf];
  vars->buf2 = new double[vars->nBuf];

  //--- Log
  char const* modeName = (vars->mode == MODE_TSR) ? "tsr  (1/B = B(1/z) / [B0(z) B0(1/z)])" :
    (vars->mode == MODE_WHITEN) ? "whiten  (1/R = 1/[B0(z) B0(1/z)], zero phase)" : "minphase  (1/B0(z), causal)";
  writer->line( "  Mode:                  %s", modeName );
  writer->line( "  Method (1/B0):         %s", vars->method == METHOD_RECURSIVE ? "recursive (B0 has same length as wavelet)" : "FIR (truncated Wiener inverse)" );
  writer->line( "  Wavelet:               %d samples (%g ms), zero time at sample %d", nb, (nb-1)*dt, shift );
  writer->line( "  White noise:           %g %%  (r0 = %g -> %g)", whiteNoise, r0orig, r[0] );
  writer->line( "  Wiener inverse length: %d samples (%g ms)", nInv, (nInv-1)*dt );
  writer->line( "  Factorization error:   %.2e   tail energy 1/B0: %.2e   recursion vs FIR: %.2e", errR, tailRatio, errH );
  writer->line( "  Wavelet phase:         %s", diffMin < 1.0e-2 ? "minimum phase (B0 = B): TSR reduces to causal 1/B" : "NOT minimum phase (B0 != B): anticausal part required" );
  writer->line( "  Minimum-phase factor B0 (first coefficients):" );
  for( int k = 0; k < nb && k < 10; k++ ) {
    writer->line( "     k=%3d   b=%12.6f   b0=%12.6f", k, wav[k], vars->b0[k] );
  }
}

//*************************************************************************************************
// Exec phase
//
void exec_mod_rm_tsr_(
  csTraceGather* traceGather,
  int* port,
  int* numTrcToKeep,
  csExecPhaseEnv* env,
  csLogWriter* writer )
{
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  csSuperHeader const* shdr = env->superHeader;
  float* samples = traceGather->trace(0)->getTraceSamples();
  int nSamples = shdr->numSamples;
  int nPad = vars->nPad;
  int nBuf = vars->nBuf;
  int nb   = vars->nb;
  int nInv = vars->nInv;
  double const* b   = vars->b;
  double const* b0  = vars->b0;
  double const* inv = vars->inv;
  double* x = vars->buf1;
  double* y = vars->buf2;

  vars->numTraces++;

  for( int i = 0; i < nBuf; i++ ) x[i] = 0.0;
  for( int i = 0; i < nSamples; i++ ) x[nPad+i] = samples[i];

  //--- 1) Correlation with b:  v_t = sum_k b_k x_{t+k}   (mode tsr only)
  if( vars->mode == MODE_TSR ) {
    for( int t = 0; t < nBuf; t++ ) {
      double sum = 0.0;
      int kmax = std::min( nb, nBuf - t );
      for( int k = 0; k < kmax; k++ ) sum += b[k] * x[t+k];
      y[t] = sum;
    }
    std::swap( x, y );
  }

  //--- 2) Causal 1/B0(z):  t = 0 -> N
  if( vars->method == METHOD_RECURSIVE ) {
    for( int t = 0; t < nBuf; t++ ) {
      double sum = x[t];
      int kmax = std::min( nb-1, t );
      for( int k = 1; k <= kmax; k++ ) sum -= b0[k] * y[t-k];
      y[t] = sum / b0[0];
    }
  }
  else {
    for( int t = 0; t < nBuf; t++ ) {
      double sum = 0.0;
      int kmax = std::min( nInv-1, t );
      for( int k = 0; k <= kmax; k++ ) sum += inv[k] * x[t-k];
      y[t] = sum;
    }
  }
  std::swap( x, y );

  //--- 3) Anticausal 1/B0(1/z):  t = N -> 0   (modes tsr, whiten)
  if( vars->mode != MODE_MINPHASE ) {
    if( vars->method == METHOD_RECURSIVE ) {
      for( int t = nBuf-1; t >= 0; t-- ) {
        double sum = x[t];
        int kmax = std::min( nb-1, nBuf-1-t );
        for( int k = 1; k <= kmax; k++ ) sum -= b0[k] * y[t+k];
        y[t] = sum / b0[0];
      }
    }
    else {
      for( int t = nBuf-1; t >= 0; t-- ) {
        double sum = 0.0;
        int kmax = std::min( nInv-1, nBuf-1-t );
        for( int k = 0; k <= kmax; k++ ) sum += inv[k] * x[t+k];
        y[t] = sum;
      }
    }
    std::swap( x, y );
  }

  //--- Output, shifted so that the wavelet zero time maps to the event time
  int shift = vars->shift;
  for( int i = 0; i < nSamples; i++ ) samples[i] = (float)x[nPad + i - shift];
}

//*************************************************************************************************
// Parameter definition (help: seaseis -m rm_tsr)
//
void params_mod_rm_tsr_( csParamDef* pdef ) {
  pdef->setModule( "RM_TSR", "Two-Sided Recursive (TSR) inverse filter of a known wavelet (RM = Roland Marschall)",
                   "Stable inverse 1/B of a wavelet that need not be minimum phase: TSR = B(1/z) / [B0(z) B0(1/z)]" );

  pdef->addDoc("The prefix RM in the module name stands for Roland Marschall.");
  pdef->addDoc("Convention: z = unit delay, B(z) = sum b_k z^k. The wavelet must be sampled at the trace sample interval.");
  pdef->addDoc("1/B(z) is causal and stable only if B is minimum phase (all roots outside the unit circle). Multiplying numerator and denominator by B(1/z) gives");
  pdef->addDoc("TSR = B(1/z) / R(z), R = B(z)B(1/z) = B0(z)B0(1/z): a correlation with b (carries the phase) followed by a causal and an anticausal recursion with the minimum-phase factor B0.");
  pdef->addDoc("B0 is obtained from the autocorrelation R by Wiener-Levinson (spike filter). The result is the unique stable inverse of B: amplitude 1/|B|, phase -phi(B).");
  pdef->addDoc("If B is already minimum phase, B0 = B and TSR reduces to the causal recursion 1/B.");
  pdef->addDoc("The output is non-causal: energy can appear before an event. Traces are zero-padded internally on both sides.");

  pdef->addParam( "wavelet", "Wavelet samples b_0 b_1 ... b_n (sampled at the trace sample interval)", NUM_VALUES_VARIABLE,
                  "Can span several lines. Use either 'wavelet' or 'filename'" );
  pdef->addValue( "", VALTYPE_NUMBER, "Wavelet sample values" );

  pdef->addParam( "filename", "ASCII file with the wavelet samples", NUM_VALUES_FIXED,
                  "Numbers separated by blanks, tabs, commas or new lines. Text after '#' is ignored" );
  pdef->addValue( "", VALTYPE_STRING, "File name" );

  pdef->addParam( "zero_time", "Time of the wavelet sample that corresponds to the event time [ms] (mode tsr)", NUM_VALUES_FIXED,
                  "0: first wavelet sample. For a zero-phase wavelet (e.g. Ricker) use the time of its peak" );
  pdef->addValue( "0", VALTYPE_NUMBER, "Zero time [ms], relative to the first wavelet sample" );

  pdef->addParam( "mode", "Filter to apply", NUM_VALUES_FIXED );
  pdef->addValue( "tsr", VALTYPE_OPTION );
  pdef->addOption( "tsr", "Full inverse 1/B = B(1/z)/[B0(z)B0(1/z)]: correct amplitude and phase" );
  pdef->addOption( "whiten", "Zero-phase two-sided filter 1/R = 1/|B|^2 only (no correlation step)" );
  pdef->addOption( "minphase", "Causal 1/B0 only (classic spike deconvolution; wrong phase if the wavelet is not minimum phase)" );

  pdef->addParam( "white_noise", "White noise added to r_0 [%]", NUM_VALUES_FIXED,
                  "Regularizes spectral zeros (roots on or close to the unit circle)" );
  pdef->addValue( "0.1", VALTYPE_NUMBER, "White noise [%]" );

  pdef->addParam( "filter_length", "Length of the Wiener inverse 1/B0 [ms]", NUM_VALUES_FIXED,
                  "Default: 4 x wavelet length. Must cover the decay of 1/B0; increase if the log reports tail energy or factorization warnings" );
  pdef->addValue( "", VALTYPE_NUMBER, "Filter length [ms]" );

  pdef->addParam( "method", "How 1/B0 is applied", NUM_VALUES_FIXED );
  pdef->addValue( "recursive", VALTYPE_OPTION );
  pdef->addOption( "recursive", "Recursion with B0 (same length as the wavelet). Falls back to 'fir' automatically if the check fails" );
  pdef->addOption( "fir", "Convolution with the truncated Wiener inverse (always stable, slower)" );
}

//************************************************************************************************
// Start exec phase
//
bool start_exec_mod_rm_tsr_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return true;
}

//************************************************************************************************
// Cleanup phase
//
void cleanup_mod_rm_tsr_( csExecPhaseEnv* env, csLogWriter* writer ) {
  VariableStruct* vars = reinterpret_cast<VariableStruct*>( env->execPhaseDef->variables() );
  writer->line( "RM_TSR: %d traces processed", vars->numTraces );
  if( vars->b    != NULL ) delete [] vars->b;
  if( vars->b0   != NULL ) delete [] vars->b0;
  if( vars->inv  != NULL ) delete [] vars->inv;
  if( vars->buf1 != NULL ) delete [] vars->buf1;
  if( vars->buf2 != NULL ) delete [] vars->buf2;
  delete vars; vars = NULL;
}

extern "C" void _params_mod_rm_tsr_( csParamDef* pdef ) {
  params_mod_rm_tsr_( pdef );
}
extern "C" void _init_mod_rm_tsr_( csParamManager* param, csInitPhaseEnv* env, csLogWriter* writer ) {
  init_mod_rm_tsr_( param, env, writer );
}
extern "C" bool _start_exec_mod_rm_tsr_( csExecPhaseEnv* env, csLogWriter* writer ) {
  return start_exec_mod_rm_tsr_( env, writer );
}
extern "C" void _exec_mod_rm_tsr_( csTraceGather* traceGather, int* port, int* numTrcToKeep, csExecPhaseEnv* env, csLogWriter* writer ) {
  exec_mod_rm_tsr_( traceGather, port, numTrcToKeep, env, writer );
}
extern "C" void _cleanup_mod_rm_tsr_( csExecPhaseEnv* env, csLogWriter* writer ) {
  cleanup_mod_rm_tsr_( env, writer );
}
