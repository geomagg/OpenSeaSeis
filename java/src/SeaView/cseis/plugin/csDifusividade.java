/* Diffuse-field tests of ambient noise along a line of receivers:
   SPAC / Bessel J0, f-k, hour-by-hour symmetry and neighbour delays (clock check)
   M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import java.util.Arrays;

/**
 * Cross-spectral matrix of the raw noise, accumulated per time block (e.g. one hour), and four tests on it.
 * <pre>
 *   S_ij(f) = sum over windows of conj(X_i) X_j          (same pre-processing as Interferometria)
 *   gamma_ij(f) = S_ij / sqrt(S_ii S_jj)                 coherency; delay of j relative to i -> phase -2 pi f tau
 *
 *   1. SPAC (Aki 1957): diffuse field  -> Re gamma(r) = J0(2 pi f r / c);  plane wave along the line -> cos(2 pi f r / c)
 *      c(f) by least squares over the distance-binned mean; also c from the first zero (kr = 2.405)
 *   2. f-k: P(k,f) = sum_ij M_ij exp(i 2 pi k (s_j - s_i)),  M = S / mean(S_ii);  k &gt; 0: energy moving towards +s
 *   3. per block: D = (P+ - P-)/(P+ + P-), P+- = energy of the f-k in the band f1-f2 and f/cmax &lt;= |k| &lt;= f/cmin;
 *      rms misfit to J0 and to cos (with c(f) of the whole record)
 *   4. per block: delay between neighbours tau = -d(phase)/d(omega) in fd1-fd2, minus the propagation part
 *      (median of tau/ds times ds). A clock or time-alignment error is constant in frequency and in time.
 * </pre>
 */
public final class csDifusividade {

  public static final class Params {
    // pre-processing (as in Interferometria)
    public double winLen = 120.0;        // [s] stacking window
    public double fmin = 0.05, fmax = 2.5;
    public int    tempNorm = 2;          // 0 none, 1 one-bit, 2 running absolute mean
    public double ramWin = 30.0;
    public boolean whiten = true;
    // blocks
    public double blockLen = 3600.0;     // [s]
    // analysis band and velocities for symmetry and delays (Scholte of line 1281: 0.15-0.45 Hz, 450-1600 m/s)
    public double f1 = 0.15, f2 = 0.45;
    public double cmin = 450, cmax = 1600;
    // band of the neighbour delays (clock check): where the coherency of close receivers is high
    public double fd1 = 0.10, fd2 = 0.35;
    // SPAC fit range
    public double cFitMin = 150, cFitMax = 4000;
    public boolean freeAmp = false;       // fit a*J0 + b instead of J0
    public int maxFreqs = 200;           // frequency groups (neighbouring FFT bins are averaged)
  }

  public static final class Result {
    public Params p;
    public int n;                        // traces
    public String[] labels;
    public double[] s;                   // position along the line [m], from the first trace in s order
    public double[][] dist;              // distance between traces [m]
    public double azimuth;               // azimuth of +s [deg from north]
    public double dx;                    // median neighbour spacing [m]
    public int[] order;                  // traces sorted by s
    public double[] freqs;               // centre of the frequency groups
    public int nblk;
    public double[] blkStart;            // [s] from the start of the traces
    public int[] blkWins;
    public int nwinTotal, nWinSamples;
    public double dfBin;
    // spectra
    public double[][] autoT;             // [g][i]
    public double[][] crossT;            // [g][2*pair]  (re, im), pair i<j
    public float[][][] autoB;            // [b][g][i]
    public float[][][] crossB;           // [b][g][2*pair]
    // 1. SPAC
    public double[] cFit, rmsJ0, rmsCos, cZero, ampA, ampB;
    public double[] cGrid;
    public double[][] spacMisfit;        // [g][c] rms misfit to J0 (NaN below the spacing limit)
    // 2. f-k
    public double[] ks;
    public double[][] fk;                // [k][g], normalized per frequency
    // 3. per block
    public double dTotal;
    public double[] dBlk, rmsJ0Blk, rmsCosBlk;
    // 4. neighbour delays (pairs order[m], order[m+1])
    public double[] tauT, excessT, ds;   // [m] in seconds / metres
    public double[][] excessB;           // [b][m]
    public double slownessT;             // median tau/ds [s/m]
    public String warnings = "";
  }

  static int pairIndex( int i, int j, int n ) {         // i < j
    return i * ( 2 * n - i - 1 ) / 2 + ( j - i - 1 );
  }

  //====================================================================
  public static Result compute( float[][] samples, double[] x, double[] y, String[] labels, double dt, Params p ) {
    int n = samples.length;
    if( n < 3 ) throw new IllegalArgumentException( "São necessários pelo menos 3 traços (receptores)" );
    int ns = samples[0].length;
    for( float[] t : samples ) ns = Math.min( ns, t.length );
    int nw = (int)Math.round( p.winLen / dt );
    if( nw < 16 ) throw new IllegalArgumentException( "Janela muito curta" );
    if( nw > ns ) throw new IllegalArgumentException( String.format( java.util.Locale.US, "A janela (%.4g s) é maior que o traço (%.4g s)", p.winLen, ns * dt ) );
    if( p.fmax > 0.5 / dt ) p.fmax = 0.5 / dt;
    if( p.f2 > p.fmax ) p.f2 = p.fmax;
    if( p.f1 < p.fmin ) p.f1 = p.fmin;
    if( p.fd2 > p.fmax ) p.fd2 = p.fmax;
    if( p.fd1 < p.fmin ) p.fd1 = p.fmin;
    if( p.fd2 <= p.fd1 ) throw new IllegalArgumentException( "Banda dos atrasos: f1 < f2, dentro da banda de pré-processamento" );
    if( p.f2 <= p.f1 ) throw new IllegalArgumentException( "Banda de análise: f1 < f2, dentro da banda de pré-processamento" );

    Result r = new Result();
    r.p = p; r.n = n; r.labels = labels;
    geometry( r, x, y );

    int nfft = 1;
    while( nfft < nw ) nfft *= 2;
    double df = 1.0 / ( nfft * dt );
    r.dfBin = df;
    int k0 = Math.max( 1, (int)Math.ceil( p.fmin / df ) ), k1 = Math.min( nfft / 2, (int)Math.floor( p.fmax / df ) );
    if( k1 < k0 + 2 ) throw new IllegalArgumentException( "Banda muito estreita para a janela escolhida" );
    int nbin = k1 - k0 + 1;
    int step = Math.max( 1, (int)Math.ceil( nbin / (double)p.maxFreqs ) );
    int nfo = ( nbin + step - 1 ) / step;
    r.freqs = new double[nfo];
    int[] cntf = new int[nfo];
    for( int k = k0; k <= k1; k++ ) { int g = ( k - k0 ) / step; r.freqs[g] += k * df; cntf[g]++; }
    for( int g = 0; g < nfo; g++ ) r.freqs[g] /= cntf[g];

    int nwin = ns / nw;
    r.nwinTotal = nwin; r.nWinSamples = nw;
    double blk = p.blockLen > 0 ? p.blockLen : ns * dt;
    r.nblk = Math.max( 1, (int)Math.ceil( nwin * nw * dt / blk - 1e-9 ) );
    int npair = n * ( n - 1 ) / 2;
    double mem = (double)r.nblk * nfo * ( n + 2.0 * npair ) * 4 + 2.0 * nfo * ( n + 2.0 * npair ) * 8;
    if( mem > 0.6 * Runtime.getRuntime().maxMemory() ) {
      throw new IllegalArgumentException( String.format( java.util.Locale.US,
          "Memória insuficiente (%.0f MB para %d traços, %d blocos, %d frequências). Aumente o bloco, reduza a banda ou use menos traços.", mem / 1e6, n, r.nblk, nfo ) );
    }
    r.autoB = new float[r.nblk][nfo][n];
    r.crossB = new float[r.nblk][nfo][2 * npair];
    r.blkStart = new double[r.nblk];
    r.blkWins = new int[r.nblk];
    for( int b = 0; b < r.nblk; b++ ) r.blkStart[b] = b * blk;

    csPluginInterferometria.Params pp = new csPluginInterferometria.Params();
    pp.fmin = p.fmin; pp.fmax = p.fmax; pp.tempNorm = p.tempNorm; pp.ramWin = p.ramWin; pp.whiten = p.whiten;
    double[] mask = csPluginInterferometria.bandMask( nfft, df, p.fmin, p.fmax );
    double[] taper = csPluginInterferometria.tukey( nw, 0.05 );
    int nram = Math.max( 1, (int)Math.round( p.ramWin / dt ) );
    double[] re = new double[nfft], im = new double[nfft];
    double[][] xr = new double[n][nbin], xi = new double[n][nbin];
    double[][] wa = new double[nfo][n], wc = new double[nfo][2 * npair];
    for( int w = 0; w < nwin; w++ ) {
      int b = Math.min( r.nblk - 1, (int)( w * nw * dt / blk ) );
      r.blkWins[b]++;
      for( int i = 0; i < n; i++ ) {
        csPluginInterferometria.preprocess( samples[i], w * nw, nw, taper, mask, nram, pp, re, im );
        for( int k = 0; k < nbin; k++ ) { xr[i][k] = re[k0 + k]; xi[i][k] = im[k0 + k]; }
      }
      for( double[] a : wa ) Arrays.fill( a, 0.0 );
      for( double[] a : wc ) Arrays.fill( a, 0.0 );
      for( int k = 0; k < nbin; k++ ) {
        int g = k / step;
        double[] ag = wa[g], cg = wc[g];
        for( int i = 0; i < n; i++ ) {
          double ar = xr[i][k], ai = xi[i][k];
          ag[i] += ar * ar + ai * ai;
          int pi = pairIndex( i, i + 1, n );
          for( int j = i + 1; j < n; j++, pi++ ) {
            double br = xr[j][k], bi = xi[j][k];
            cg[2 * pi]     += ar * br + ai * bi;        // conj(Xi) Xj
            cg[2 * pi + 1] += ar * bi - ai * br;
          }
        }
      }
      for( int g = 0; g < nfo; g++ ) {
        float[] ab = r.autoB[b][g], cb = r.crossB[b][g];
        for( int i = 0; i < n; i++ ) ab[i] += (float)wa[g][i];
        for( int q = 0; q < 2 * npair; q++ ) cb[q] += (float)wc[g][q];
      }
    }
    r.autoT = new double[nfo][n];
    r.crossT = new double[nfo][2 * npair];
    for( int b = 0; b < r.nblk; b++ ) for( int g = 0; g < nfo; g++ ) {
      for( int i = 0; i < n; i++ ) r.autoT[g][i] += r.autoB[b][g][i];
      for( int q = 0; q < 2 * npair; q++ ) r.crossT[g][q] += r.crossB[b][g][q];
    }
    if( nwin * nw < ns ) r.warnings += String.format( java.util.Locale.US, "Últimos %.0f s não usados (janela incompleta). ", ( ns - nwin * nw ) * dt );
    if( r.nblk > 1 && r.blkWins[r.nblk - 1] < 0.5 * r.blkWins[0] ) r.warnings += "O último bloco é mais curto que os outros. ";

    spac( r );
    fk( r );
    blocks( r );
    delays( r );
    return r;
  }

  //--------------------------------------------------------------------
  static void geometry( Result r, double[] x, double[] y ) {
    int n = r.n;
    double mx = 0, my = 0;
    for( int i = 0; i < n; i++ ) { mx += x[i]; my += y[i]; }
    mx /= n; my /= n;
    double sxx = 0, syy = 0, sxy = 0;
    for( int i = 0; i < n; i++ ) { double a = x[i] - mx, b = y[i] - my; sxx += a * a; syy += b * b; sxy += a * b; }
    double th = 0.5 * Math.atan2( 2 * sxy, sxx - syy );
    double ux = Math.cos( th ), uy = Math.sin( th );
    // +s from the first to the last trace of the pane
    if( ( x[n - 1] - x[0] ) * ux + ( y[n - 1] - y[0] ) * uy < 0 ) { ux = -ux; uy = -uy; }
    r.s = new double[n];
    double smin = Double.POSITIVE_INFINITY;
    for( int i = 0; i < n; i++ ) { r.s[i] = ( x[i] - mx ) * ux + ( y[i] - my ) * uy; smin = Math.min( smin, r.s[i] ); }
    for( int i = 0; i < n; i++ ) r.s[i] -= smin;
    r.azimuth = ( Math.toDegrees( Math.atan2( ux, uy ) ) + 360.0 ) % 360.0;
    r.dist = new double[n][n];
    for( int i = 0; i < n; i++ ) for( int j = 0; j < n; j++ ) r.dist[i][j] = Math.hypot( x[i] - x[j], y[i] - y[j] );
    Integer[] o = new Integer[n];
    for( int i = 0; i < n; i++ ) o[i] = i;
    Arrays.sort( o, ( a, b ) -> Double.compare( r.s[a], r.s[b] ) );
    r.order = new int[n];
    for( int i = 0; i < n; i++ ) r.order[i] = o[i];
    double[] d = new double[n - 1];
    for( int m = 0; m < n - 1; m++ ) d[m] = r.dist[r.order[m]][r.order[m + 1]];
    Arrays.sort( d );
    r.dx = d[( n - 1 ) / 2];
    if( !( r.dx > 0 ) ) throw new IllegalArgumentException( "Receptores com a mesma posição (rec_x, rec_y): confira as coordenadas" );
  }

  /** Coherency of pair (i,j) at group g; out = {re, im} */
  static void coh( double[] auto, double[] cross, int i, int j, int n, double[] out ) {
    double a = Math.sqrt( auto[i] * auto[j] );
    if( !( a > 0 ) ) { out[0] = Double.NaN; out[1] = Double.NaN; return; }
    int a0 = i, b0 = j; double sg = 1;
    if( i > j ) { a0 = j; b0 = i; sg = -1; }            // gamma_ji = conj(gamma_ij)
    int q = pairIndex( a0, b0, n );
    out[0] = cross[2 * q] / a; out[1] = sg * cross[2 * q + 1] / a;
  }

  /** Mean Re gamma (and Im) per distance bin of width dx. Returns {rMean[], reMean[], imMean[], count[]} */
  static double[][] binned( Result r, double[] auto, double[] cross ) {
    int n = r.n;
    double rmax = 0;
    for( int i = 0; i < n; i++ ) for( int j = i + 1; j < n; j++ ) rmax = Math.max( rmax, r.dist[i][j] );
    int nb = (int)Math.round( rmax / r.dx ) + 2;
    double[] rs = new double[nb], res = new double[nb], ims = new double[nb], cnt = new double[nb];
    double[] c = new double[2];
    for( int i = 0; i < n; i++ ) for( int j = i + 1; j < n; j++ ) {
      coh( auto, cross, i, j, n, c );
      if( Double.isNaN( c[0] ) ) continue;
      int b = (int)Math.round( r.dist[i][j] / r.dx );
      // Im with the sign of s_j - s_i: wave moving towards +s gives Im = -sin(kr) (< 0 for small kr)
      double sg = Math.signum( r.s[j] - r.s[i] );
      rs[b] += r.dist[i][j]; res[b] += c[0]; ims[b] += c[1] * ( sg == 0 ? 1 : sg ); cnt[b]++;
    }
    int m = 0;
    for( int b = 0; b < nb; b++ ) if( cnt[b] > 0 ) m++;
    double[][] o = new double[4][m];
    m = 0;
    for( int b = 0; b < nb; b++ ) if( cnt[b] > 0 ) { o[0][m] = rs[b] / cnt[b]; o[1][m] = res[b] / cnt[b]; o[2][m] = ims[b] / cnt[b]; o[3][m] = cnt[b]; m++; }
    return o;
  }

  /**
   * rms of (binned mean - model), weighted by the number of pairs; model 0 = J0, 1 = cos.
   * free = true: model a*J0(kr) + b with a, b by least squares (0 &lt; a &lt;= 1.2); b absorbs energy arriving
   * vertically (k ~ 0, adds a positive constant) and a the incoherent noise (lowers the amplitude).
   * out (optional) receives {a, b}
   */
  static double misfit( double[][] bins, double f, double c, int model, boolean free, double[] out ) {
    int nb = bins[0].length;
    double sw = 0, sm = 0, sy = 0, smm = 0, smy = 0;
    double[] m = new double[nb];
    for( int b = 0; b < nb; b++ ) {
      double kr = 2 * Math.PI * f * bins[0][b] / c;
      m[b] = model == 0 ? j0( kr ) : Math.cos( kr );
      double w = bins[3][b], yv = bins[1][b];
      sw += w; sm += w * m[b]; sy += w * yv; smm += w * m[b] * m[b]; smy += w * m[b] * yv;
    }
    if( !( sw > 0 ) ) return Double.NaN;
    double a = 1, bb = 0;
    if( free ) {
      double den = sw * smm - sm * sm;
      if( Math.abs( den ) > 1e-12 * sw * sw ) {
        a = ( sw * smy - sm * sy ) / den;
        a = Math.max( 0.05, Math.min( 1.2, a ) );
      }
      bb = ( sy - a * sm ) / sw;
    }
    double s = 0;
    for( int b = 0; b < nb; b++ ) { double e = bins[1][b] - ( a * m[b] + bb ); s += bins[3][b] * e * e; }
    if( out != null ) { out[0] = a; out[1] = bb; }
    return Math.sqrt( s / sw );
  }
  static double misfit( double[][] bins, double f, double c, int model ) { return misfit( bins, f, c, model, false, null ); }

  static double[] cGrid( Params p ) {
    int nc = 300;
    double[] c = new double[nc];
    for( int k = 0; k < nc; k++ ) c[k] = p.cFitMin * Math.pow( p.cFitMax / p.cFitMin, k / (double)( nc - 1 ) );
    return c;
  }
  /** lowest c whose first J0 zero is beyond one receiver spacing (below it J0 ~ 0 at every measured distance) */
  static double cLow( double f, double dx ) { return 2 * Math.PI * f * dx / 2.404826; }

  /**
   * Best c on the grid for J0; returns {c, rmsJ0, rmsCos(same c), a, b}; misfitOut (optional) receives the
   * J0 misfit for every c of the grid (NaN below cLow)
   */
  static double[] fitC( double[][] bins, double f, double[] cg, double dx, boolean free, double[] misfitOut ) {
    double[] none = { Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN };
    double lo = cLow( f, dx );
    double best = Double.POSITIVE_INFINITY;
    int kb = -1, k0 = -1;
    for( int k = 0; k < cg.length; k++ ) {
      if( cg[k] < lo ) { if( misfitOut != null ) misfitOut[k] = Double.NaN; continue; }
      if( k0 < 0 ) k0 = k;
      double e = misfit( bins, f, cg[k], 0, free, null );
      if( misfitOut != null ) misfitOut[k] = e;
      if( e < best ) { best = e; kb = k; }
    }
    if( kb < 0 ) return none;
    if( kb == k0 || kb == cg.length - 1 ) { none[1] = best; return none; }   // minimum at a bound: no fit
    double[] ab = new double[2];
    misfit( bins, f, cg[kb], 0, free, ab );
    return new double[]{ cg[kb], best, misfit( bins, f, cg[kb], 1, free, null ), ab[0], ab[1] };
  }

  /** c from the first + to - crossing of the binned mean: kr0 = 2.405 */
  static double zeroC( double[][] bins, double f ) {
    for( int b = 1; b < bins[0].length; b++ ) {
      double a = bins[1][b - 1], c = bins[1][b];
      if( a > 0 && c <= 0 ) {
        double r0 = bins[0][b - 1] + ( bins[0][b] - bins[0][b - 1] ) * a / ( a - c );
        return 2 * Math.PI * f * r0 / 2.404826;
      }
    }
    return Double.NaN;
  }

  static void spac( Result r ) {
    int nf = r.freqs.length;
    r.cFit = new double[nf]; r.rmsJ0 = new double[nf]; r.rmsCos = new double[nf]; r.cZero = new double[nf];
    r.ampA = new double[nf]; r.ampB = new double[nf];
    r.cGrid = cGrid( r.p );
    r.spacMisfit = new double[nf][r.cGrid.length];
    for( int g = 0; g < nf; g++ ) {
      double[][] bins = binned( r, r.autoT[g], r.crossT[g] );
      double[] fc = fitC( bins, r.freqs[g], r.cGrid, r.dx, r.p.freeAmp, r.spacMisfit[g] );
      r.cFit[g] = fc[0]; r.rmsJ0[g] = fc[1]; r.rmsCos[g] = fc[2]; r.ampA[g] = fc[3]; r.ampB[g] = fc[4];
      r.cZero[g] = zeroC( bins, r.freqs[g] );
    }
  }

  //--------------------------------------------------------------------
  static double kNyq( Result r ) { return 0.5 / r.dx; }

  /** f-k power of the frequency group g, for the wavenumbers ks (not normalized; P >= 0) */
  static double[] beam( Result r, double[] auto, double[] cross, double[] ks, double[][] cs, double[][] sn ) {
    int n = r.n;
    double[] P = new double[ks.length];
    double m = 0;
    for( int i = 0; i < n; i++ ) m += auto[i];
    m /= n;
    if( !( m > 0 ) ) return P;
    double d0 = n;                                       // sum of auto / m
    for( int ik = 0; ik < ks.length; ik++ ) {
      double v = d0;
      int q = 0;
      for( int i = 0; i < n; i++ ) for( int j = i + 1; j < n; j++, q++ ) {
        // 2 Re( S_ij exp(i 2 pi k (s_j - s_i)) )
        v += 2 * ( cross[2 * q] * cs[ik][q] - cross[2 * q + 1] * sn[ik][q] ) / m;
      }
      P[ik] = Math.max( 0.0, v );
    }
    return P;
  }

  static double[][][] phaseTables( Result r, double[] ks ) {
    int n = r.n, npair = n * ( n - 1 ) / 2;
    double[][] cs = new double[ks.length][npair], sn = new double[ks.length][npair];
    for( int ik = 0; ik < ks.length; ik++ ) {
      int q = 0;
      for( int i = 0; i < n; i++ ) for( int j = i + 1; j < n; j++, q++ ) {
        double a = 2 * Math.PI * ks[ik] * ( r.s[j] - r.s[i] );
        cs[ik][q] = Math.cos( a ); sn[ik][q] = Math.sin( a );
      }
    }
    return new double[][][]{ cs, sn };
  }

  static void fk( Result r ) {
    int nk = 201;
    double kmax = kNyq( r );
    r.ks = new double[nk];
    for( int k = 0; k < nk; k++ ) r.ks[k] = -kmax + 2 * kmax * k / ( nk - 1 );
    double[][][] t = phaseTables( r, r.ks );
    int nf = r.freqs.length;
    r.fk = new double[nk][nf];
    for( int g = 0; g < nf; g++ ) {
      double[] P = beam( r, r.autoT[g], r.crossT[g], r.ks, t[0], t[1] );
      double mx = 0;
      for( double v : P ) mx = Math.max( mx, v );
      for( int k = 0; k < nk; k++ ) r.fk[k][g] = mx > 0 ? P[k] / mx : 0;
    }
  }

  /** D = (P+ - P-)/(P+ + P-) in f1-f2 and f/cmax <= |k| <= f/cmin (wavenumbers on a fine grid up to Nyquist) */
  static double symmetry( Result r, double[] auto0, double[][] autoG, double[][] crossG, double[] ks, double[][][] t ) {
    double pp = 0, pm = 0;
    for( int g = 0; g < r.freqs.length; g++ ) {
      double f = r.freqs[g];
      if( f < r.p.f1 || f > r.p.f2 ) continue;
      double[] P = beam( r, autoG[g], crossG[g], ks, t[0], t[1] );
      double tot = 0, a = 0, b = 0;
      for( int k = 0; k < ks.length; k++ ) {
        double ak = Math.abs( ks[k] );
        tot += P[k];
        if( ak < f / r.p.cmax || ak > f / r.p.cmin ) continue;
        if( ks[k] > 0 ) a += P[k]; else b += P[k];
      }
      if( tot > 0 ) { pp += a / tot; pm += b / tot; }   // each frequency with the same weight
    }
    return ( pp + pm > 0 ) ? ( pp - pm ) / ( pp + pm ) : Double.NaN;
  }

  static void blocks( Result r ) {
    int nk = 161;
    double kmax = kNyq( r );
    double[] ks = new double[nk];
    for( int k = 0; k < nk; k++ ) ks[k] = -kmax + 2 * kmax * k / ( nk - 1 );
    double[][][] t = phaseTables( r, ks );
    r.dTotal = symmetry( r, null, r.autoT, r.crossT, ks, t );
    r.dBlk = new double[r.nblk]; r.rmsJ0Blk = new double[r.nblk]; r.rmsCosBlk = new double[r.nblk];
    int nf = r.freqs.length;
    for( int b = 0; b < r.nblk; b++ ) {
      double[][] ag = new double[nf][], cg = new double[nf][];
      for( int g = 0; g < nf; g++ ) { ag[g] = toD( r.autoB[b][g] ); cg[g] = toD( r.crossB[b][g] ); }
      r.dBlk[b] = r.blkWins[b] > 0 ? symmetry( r, null, ag, cg, ks, t ) : Double.NaN;
      double sj = 0, sc = 0; int m = 0;
      for( int g = 0; g < nf; g++ ) {
        double f = r.freqs[g];
        if( f < r.p.f1 || f > r.p.f2 || Double.isNaN( r.cFit[g] ) || r.blkWins[b] == 0 ) continue;
        double[][] bins = binned( r, ag[g], cg[g] );
        double ej = misfit( bins, f, r.cFit[g], 0, r.p.freeAmp, null );
        double ec = misfit( bins, f, r.cFit[g], 1, r.p.freeAmp, null );
        if( Double.isNaN( ej ) ) continue;
        sj += ej; sc += ec; m++;
      }
      r.rmsJ0Blk[b] = m > 0 ? sj / m : Double.NaN;
      r.rmsCosBlk[b] = m > 0 ? sc / m : Double.NaN;
    }
  }

  static double[] toD( float[] a ) {
    double[] d = new double[a.length];
    for( int i = 0; i < a.length; i++ ) d[i] = a[i];
    return d;
  }

  //--------------------------------------------------------------------
  /** Delay of j relative to i [s] from the slope of the coherency phase in f1-f2 (weights |gamma|^2) */
  static double delay( Result r, double[][] autoG, double[][] crossG, int i, int j ) {
    double sw = 0, sx = 0, sy = 0, sxx = 0, sxy = 0, prev = Double.NaN, off = 0;
    double[] c = new double[2];
    for( int g = 0; g < r.freqs.length; g++ ) {
      double f = r.freqs[g];
      if( f < r.p.fd1 || f > r.p.fd2 ) continue;
      coh( autoG[g], crossG[g], i, j, r.n, c );
      if( Double.isNaN( c[0] ) ) continue;
      double ph = Math.atan2( c[1], c[0] );
      if( !Double.isNaN( prev ) ) {                       // unwrap
        while( ph + off - prev > Math.PI ) off -= 2 * Math.PI;
        while( ph + off - prev < -Math.PI ) off += 2 * Math.PI;
      }
      ph += off; prev = ph;
      double w = c[0] * c[0] + c[1] * c[1], om = 2 * Math.PI * f;
      sw += w; sx += w * om; sy += w * ph; sxx += w * om * om; sxy += w * om * ph;
    }
    double den = sw * sxx - sx * sx;
    if( !( sw > 0 ) || Math.abs( den ) < 1e-30 ) return Double.NaN;
    double slope = ( sw * sxy - sx * sy ) / den;
    return -slope;
  }

  static void delays( Result r ) {
    int m = r.n - 1;
    r.tauT = new double[m]; r.excessT = new double[m]; r.ds = new double[m];
    double[] sl = new double[m];
    for( int k = 0; k < m; k++ ) {
      int i = r.order[k], j = r.order[k + 1];
      r.ds[k] = r.s[j] - r.s[i];
      r.tauT[k] = delay( r, r.autoT, r.crossT, i, j );
      sl[k] = r.tauT[k] / r.ds[k];
    }
    r.slownessT = medianNaN( sl );
    for( int k = 0; k < m; k++ ) r.excessT[k] = r.tauT[k] - r.slownessT * r.ds[k];
    r.excessB = new double[r.nblk][m];
    int nf = r.freqs.length;
    for( int b = 0; b < r.nblk; b++ ) {
      double[][] ag = new double[nf][], cg = new double[nf][];
      for( int g = 0; g < nf; g++ ) { ag[g] = toD( r.autoB[b][g] ); cg[g] = toD( r.crossB[b][g] ); }
      double[] tb = new double[m], sb = new double[m];
      for( int k = 0; k < m; k++ ) {
        tb[k] = r.blkWins[b] > 0 ? delay( r, ag, cg, r.order[k], r.order[k + 1] ) : Double.NaN;
        sb[k] = tb[k] / r.ds[k];
      }
      double slb = medianNaN( sb );
      for( int k = 0; k < m; k++ ) r.excessB[b][k] = tb[k] - slb * r.ds[k];
    }
  }

  static double medianNaN( double[] v ) {
    double[] a = Arrays.stream( v ).filter( d -> !Double.isNaN( d ) ).sorted().toArray();
    if( a.length == 0 ) return Double.NaN;
    return ( a.length % 2 == 1 ) ? a[a.length / 2] : 0.5 * ( a[a.length / 2 - 1] + a[a.length / 2] );
  }
  /** robust spread: 1.4826 * median |x - median| */
  static double mad( double[] v ) {
    double med = medianNaN( v );
    double[] d = new double[v.length];
    for( int i = 0; i < v.length; i++ ) d[i] = Math.abs( v[i] - med );
    return 1.4826 * medianNaN( d );
  }

  /** Suspect neighbour pairs: |excess| above max(60 ms, 3 sigma) in the whole record and same sign in at least 75% of the blocks */
  static boolean[] suspects( Result r, double[] thrOut ) {
    int m = r.excessT.length;
    double sig = mad( r.excessT );
    double thr = Math.max( 0.060, 3 * sig );
    thrOut[0] = thr;
    boolean[] s = new boolean[m];
    for( int k = 0; k < m; k++ ) {
      double e = r.excessT[k];
      if( Double.isNaN( e ) || Math.abs( e ) < thr ) continue;
      int same = 0, tot = 0;
      for( int b = 0; b < r.nblk; b++ ) {
        double v = r.excessB[b][k];
        if( Double.isNaN( v ) ) continue;
        tot++;
        if( Math.signum( v ) == Math.signum( e ) && Math.abs( v ) > 0.5 * thr ) same++;
      }
      s[k] = tot == 0 || same >= 0.75 * tot;
    }
    return s;
  }

  //--------------------------------------------------------------------
  /** Bessel function of the first kind, order 0 (rational approximations, |error| &lt; 1e-7) */
  public static double j0( double x ) {
    double ax = Math.abs( x );
    if( ax < 8.0 ) {
      double y = x * x;
      double a1 = 57568490574.0 + y * ( -13362590354.0 + y * ( 651619640.7 + y * ( -11214424.18 + y * ( 77392.33017 + y * ( -184.9052456 ) ) ) ) );
      double a2 = 57568490411.0 + y * ( 1029532985.0 + y * ( 9494680.718 + y * ( 59272.64853 + y * ( 267.8532712 + y ) ) ) );
      return a1 / a2;
    }
    double z = 8.0 / ax, y = z * z, xx = ax - 0.785398164;
    double a1 = 1.0 + y * ( -0.1098628627e-2 + y * ( 0.2734510407e-4 + y * ( -0.2073370639e-5 + y * 0.2093887211e-6 ) ) );
    double a2 = -0.1562499995e-1 + y * ( 0.1430488765e-3 + y * ( -0.6911147651e-5 + y * ( 0.7621095161e-6 - y * 0.934935152e-7 ) ) );
    return Math.sqrt( 0.636619772 / ax ) * ( Math.cos( xx ) * a1 - z * Math.sin( xx ) * a2 );
  }
}
