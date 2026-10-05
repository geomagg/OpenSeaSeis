/* SeaView clock QC by noise cross-correlation: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import cseis.seis.csHeader;
import cseis.seis.csHeaderDef;
import java.time.LocalDateTime;
import java.time.ZoneOffset;
import java.util.Arrays;

/**
 * Clock QC of OBN nodes from ambient-noise cross-correlation (pure computation, no Swing).
 * <p>
 * Steps:
 * <ol>
 * <li><b>Alignment</b>: every trace has its own start time (trace headers time_samp1, or
 *     time_year/day/hour/min/sec). All traces are cut to the COMMON time window
 *     [latest start, earliest end] -- a trace starting at 7 s and one starting at 12 s both start
 *     at 12 s, and all end where the trace that started at 7 s ends. Sub-sample shifts are
 *     applied by linear interpolation.</li>
 * <li><b>Pre-processing</b>, window by window (window length W): demean, taper, band-pass
 *     (cosine-tapered mask in frequency), optional one-bit normalization and spectral whitening.</li>
 * <li><b>Cross-correlation</b> of every node with its K next neighbours (node order), stacked over
 *     all windows (and, separately, per period, to follow the clock drift).</li>
 * <li><b>Clock measurement</b> from the time symmetry of each stacked correlation C(tau): for
 *     a diffuse noise field C is symmetric; a clock error shifts it. With a = e_b - e_a (clock
 *     errors of the two nodes), the positive branch P(k) = C(k) equals the negative branch
 *     N(k) = C(-k) shifted by 2a: P(k) = N(k + 2a). The lag of the maximum of the correlation
 *     between P and N gives 2a. The amplitude ratio of the two branches is reported too.</li>
 * <li><b>Per-node clock error</b>: least squares of e_b - e_a = a over all pairs (weight = quality
 *     of the symmetry fit), with sum(e) = 0 or e(reference node) = 0; per period too, giving a
 *     linear clock drift [ms/day].</li>
 * </ol>
 * Sign: e &gt; 0 means the node's clock is LATE (behind true time): an event that happened at true
 * time T is stamped T - e, i.e. it appears e earlier in the trace than it should.
 */
public final class csClockXcorr {

  private csClockXcorr() {}

  //====================================================================
  public static final class Params {
    public double fmin = 0.2;          // [Hz]
    public double fmax = 2.0;          // [Hz]
    public double winSec = 120.0;      // correlation window [s]
    public double maxLagSec = 20.0;    // max lag kept [s]
    public double minLagSec = 0.2;     // lags |tau| < this are ignored in the symmetry fit [s]
    public double periodMin = 60.0;    // stacking period for the drift [min], <= 0: only total
    public boolean oneBit = true;
    public boolean whiten = true;
    public int neighbors = 5;          // pairs: each node with its next K nodes
    public double minQuality = 0.2;    // pairs with symmetry correlation below this are not used
    public int refNode = Integer.MIN_VALUE;  // node with e = 0; MIN_VALUE: mean(e) = 0
  }

  public interface Progress {
    /** @return false to cancel */
    boolean update( int done, int total, String message );
  }

  /** Common time window of the traces */
  public static final class Alignment {
    public double tStart;        // [s] absolute (latest start)
    public double tEnd;          // [s] absolute (earliest end)
    public int nCommon;          // samples in the common window
    public double[] shift;       // per trace: index of the common start in the trace [samples, fractional]
    public double[] t0;          // per trace: original start [s]
  }

  public static final class Result {
    public Alignment align;
    public int[] nodes;            // sorted node numbers (one trace each)
    public double dt;              // [s]
    public int nl;                 // lags -nl..nl
    public int nWindows;
    public int nPeriods;
    public double[] periodTime;    // centre of each period, seconds after tStart
    // pairs
    public int[] pairA, pairB;     // indices into nodes
    public float[][] cc;           // total stack per pair, 2*nl+1 samples (index nl = lag 0)
    public double[] pairShift;     // a = e_b - e_a [s] (NaN: not measured)
    public double[] pairQuality;   // symmetry correlation coefficient
    public double[] pairAmpPos, pairAmpNeg;   // max |C| on the positive / negative branch
    // nodes
    public double[] nodeErr;           // clock error e [s] from the total stack (NaN: no valid pair)
    public double[][] nodeErrPeriod;   // [node][period] [s]
    public double[] nodeDrift;         // [s/day] (NaN: fewer than 2 periods)
    public double[] nodeResid;         // nodeErr minus a robust linear trend along the line [s]
    public boolean[] nodeSuspect;      // |nodeResid| clearly above the scatter of the other nodes
    public int[] nodePairsUsed;
  }

  //====================================================================
  // Time headers
  //====================================================================
  public static int headerIndex( csHeaderDef[] defs, String name ) {
    if( defs == null ) return -1;
    for( int i = 0; i < defs.length; i++ ) {
      if( defs[i] != null && defs[i].name != null && defs[i].name.equalsIgnoreCase( name ) ) return i;
    }
    return -1;
  }
  private static double hv( csHeader[] h, int i ) {
    return ( i >= 0 && h != null && i < h.length && h[i] != null ) ? h[i].doubleValue() : 0.0;
  }
  /**
   * Absolute start time of a trace [s since 1970-01-01 UTC], from time_samp1 (+time_samp1_us), or
   * time_year/time_day (julian)/time_hour/time_min/time_sec (+time_msec/time_usec). Without
   * time_year, year 1970 is used (fine for differences between traces of the same record).
   * @return NaN if the trace has none of these headers
   */
  public static double startTimeSec( csHeader[] h, csHeaderDef[] defs ) {
    int iS = headerIndex( defs, "time_samp1" );
    if( iS >= 0 && hv( h, iS ) != 0.0 ) {
      return hv( h, iS ) + 1.0e-6 * hv( h, headerIndex( defs, "time_samp1_us" ) );
    }
    int iD = headerIndex( defs, "time_day" ), iH = headerIndex( defs, "time_hour" );
    int iM = headerIndex( defs, "time_min" ), iSec = headerIndex( defs, "time_sec" );
    if( iD < 0 && iH < 0 && iM < 0 && iSec < 0 ) return Double.NaN;
    int year = (int)hv( h, headerIndex( defs, "time_year" ) );
    if( year <= 0 ) year = 1970;
    else if( year < 100 ) year += 2000;
    int day = Math.max( 1, (int)hv( h, iD ) );
    long epoch = LocalDateTime.of( year, 1, 1, 0, 0 ).toEpochSecond( ZoneOffset.UTC ) + 86400L * ( day - 1 );
    return epoch + 3600.0 * hv( h, iH ) + 60.0 * hv( h, iM ) + hv( h, iSec )
         + 1.0e-3 * hv( h, headerIndex( defs, "time_msec" ) ) + 1.0e-6 * hv( h, headerIndex( defs, "time_usec" ) );
  }

  /**
   * Common window of traces with start times t0 [s] and ns samples each, sample interval dt [s].
   * @return null if the traces do not overlap in time
   */
  public static Alignment align( double[] t0, int ns, double dt ) {
    Alignment a = new Alignment();
    a.t0 = t0.clone();
    a.tStart = Double.NEGATIVE_INFINITY;
    a.tEnd = Double.POSITIVE_INFINITY;
    for( double t : t0 ) {
      a.tStart = Math.max( a.tStart, t );
      a.tEnd = Math.min( a.tEnd, t + ( ns - 1 ) * dt );
    }
    if( !( a.tEnd > a.tStart ) ) return null;
    a.nCommon = (int)Math.floor( ( a.tEnd - a.tStart ) / dt + 1.0e-6 ) + 1;
    a.shift = new double[t0.length];
    for( int i = 0; i < t0.length; i++ ) a.shift[i] = ( a.tStart - t0[i] ) / dt;
    // keep every interpolated sample inside the trace
    int maxN = Integer.MAX_VALUE;
    for( int i = 0; i < t0.length; i++ ) maxN = Math.min( maxN, (int)Math.floor( ( ns - 1 ) - a.shift[i] + 1.0e-6 ) + 1 );
    a.nCommon = Math.min( a.nCommon, maxN );
    return a;
  }

  /** Sample of trace s at fractional index x (linear interpolation, 0 outside) */
  public static double interp( float[] s, double x ) {
    if( x < 0.0 || x > s.length - 1 ) return 0.0;
    int i = (int)x;
    if( i >= s.length - 1 ) return s[s.length - 1];
    double w = x - i;
    return ( 1.0 - w ) * s[i] + w * s[i + 1];
  }

  //====================================================================
  // Main computation
  //====================================================================
  /**
   * @param nodes   node number of each trace (one trace per node, any order)
   * @param samples samples of each trace (same length)
   * @param t0      absolute start time of each trace [s]
   * @param dt      sample interval [s]
   * @return result, or null if cancelled
   */
  public static Result run( int[] nodes, float[][] samples, double[] t0, double dt, Params p, Progress progress ) {
    int n = nodes.length;
    if( n < 2 ) throw new IllegalArgumentException( "Need at least 2 nodes" );
    // sort by node number
    Integer[] order = new Integer[n];
    for( int i = 0; i < n; i++ ) order[i] = i;
    Arrays.sort( order, ( x, y ) -> Integer.compare( nodes[x], nodes[y] ) );
    int[] nd = new int[n];
    float[][] s = new float[n][];
    double[] t = new double[n];
    for( int i = 0; i < n; i++ ) {
      nd[i] = nodes[order[i]];
      s[i] = samples[order[i]];
      t[i] = t0[order[i]];
    }
    int ns = s[0].length;
    Alignment al = align( t, ns, dt );
    if( al == null ) throw new IllegalArgumentException( "The traces have no common time window (no overlap between start and end times)" );

    int nw = (int)Math.round( p.winSec / dt );
    if( nw < 16 ) throw new IllegalArgumentException( "Correlation window too short" );
    int nWin = al.nCommon / nw;
    if( nWin < 1 ) throw new IllegalArgumentException( String.format( java.util.Locale.US,
        "Common window (%.1f s) shorter than one correlation window (%.1f s)", al.nCommon * dt, p.winSec ) );
    int nl = (int)Math.round( p.maxLagSec / dt );
    nl = Math.min( nl, nw - 1 );
    int kmin = Math.max( 0, (int)Math.round( p.minLagSec / dt ) );
    if( kmin >= nl - 2 ) throw new IllegalArgumentException( "Minimum lag must be smaller than the maximum lag" );
    int nfft = 1;
    while( nfft < 2 * nw ) nfft <<= 1;

    int winsPerPeriod = ( p.periodMin > 0.0 ) ? Math.max( 1, (int)Math.round( p.periodMin * 60.0 / p.winSec ) ) : nWin;
    int nPer = (int)Math.ceil( nWin / (double)winsPerPeriod );

    // pairs
    int K = Math.max( 1, Math.min( p.neighbors, n - 1 ) );
    int nPairs = 0;
    for( int i = 0; i < n; i++ ) nPairs += Math.min( K, n - 1 - i );
    int[] pa = new int[nPairs], pb = new int[nPairs];
    int[][] pairIndex = new int[n][n];
    for( int[] row : pairIndex ) Arrays.fill( row, -1 );
    int ip = 0;
    for( int i = 0; i < n; i++ ) {
      for( int j = i + 1; j <= Math.min( n - 1, i + K ); j++ ) {
        pa[ip] = i; pb[ip] = j; pairIndex[i][j] = ip; ip++;
      }
    }

    int lenCC = 2 * nl + 1;
    float[][] total = new float[nPairs][lenCC];
    double[][] period = new double[nPairs][lenCC];
    double[][] perShift = new double[nPer][nPairs];
    double[][] perQual = new double[nPer][nPairs];

    // frequency mask (cosine taper over 20% of the band at each edge)
    double df = 1.0 / ( nfft * dt );
    double[] mask = bandMask( nfft, df, p.fmin, p.fmax );
    double[] taper = tukey( nw, 0.05 );

    double[][] specRe = new double[n][nfft], specIm = new double[n][nfft];
    double[] cr = new double[nfft], ci = new double[nfft];
    int done = 0;
    int totalSteps = nWin * n;
    int iPer = 0, winInPer = 0;
    for( int w = 0; w < nWin; w++ ) {
      // spectra of all nodes for this window (arrays reused from window to window)
      for( int i = 0; i < n; i++ ) {
        preprocess( s[i], al.shift[i] + (double)w * nw, nw, taper, mask, p, specRe[i], specIm[i] );
        done++;
        if( progress != null && ( done % 8 == 0 ) && !progress.update( done, totalSteps,
            String.format( java.util.Locale.US, "Janela %d/%d", w + 1, nWin ) ) ) return null;
      }
      for( int k = 0; k < nPairs; k++ ) {
        int a = pa[k], b = pb[k];
        // C_ab(tau) = sum u_a(t) u_b(t+tau)  <->  conj(A) * B
        for( int f = 0; f < nfft; f++ ) {
          double ar = specRe[a][f], ai = specIm[a][f], br = specRe[b][f], bi = specIm[b][f];
          cr[f] = ar * br + ai * bi;
          ci[f] = ar * bi - ai * br;
        }
        ifft( cr, ci );
        for( int l = -nl; l <= nl; l++ ) {
          double v = cr[ l >= 0 ? l : nfft + l ];
          total[k][l + nl] += (float)v;
          period[k][l + nl] += v;
        }
      }
      winInPer++;
      if( winInPer == winsPerPeriod || w == nWin - 1 ) {
        for( int k = 0; k < nPairs; k++ ) {
          double[] m = measure( period[k], nl, kmin );
          perShift[iPer][k] = m[0] * dt;
          perQual[iPer][k] = m[1];
          Arrays.fill( period[k], 0.0 );
        }
        iPer++;
        winInPer = 0;
      }
    }

    Result r = new Result();
    r.align = al;
    r.nodes = nd;
    r.dt = dt;
    r.nl = nl;
    r.nWindows = nWin;
    r.nPeriods = nPer;
    r.periodTime = new double[nPer];
    for( int q = 0; q < nPer; q++ ) {
      int w0 = q * winsPerPeriod, w1 = Math.min( nWin, w0 + winsPerPeriod );
      r.periodTime[q] = 0.5 * ( w0 + w1 ) * nw * dt;
    }
    r.pairA = pa; r.pairB = pb; r.cc = total;
    r.pairShift = new double[nPairs];
    r.pairQuality = new double[nPairs];
    r.pairAmpPos = new double[nPairs];
    r.pairAmpNeg = new double[nPairs];
    for( int k = 0; k < nPairs; k++ ) {
      double[] c = new double[lenCC];
      for( int l = 0; l < lenCC; l++ ) c[l] = total[k][l];
      double[] m = measure( c, nl, kmin );
      r.pairShift[k] = m[0] * dt;
      r.pairQuality[k] = m[1];
      r.pairAmpPos[k] = m[2];
      r.pairAmpNeg[k] = m[3];
    }
    int refIdx = -1;
    if( p.refNode != Integer.MIN_VALUE ) {
      for( int i = 0; i < n; i++ ) if( nd[i] == p.refNode ) refIdx = i;
    }
    r.nodePairsUsed = new int[n];
    r.nodeErr = solve( n, pa, pb, r.pairShift, r.pairQuality, p.minQuality, refIdx, r.nodePairsUsed );
    r.nodeErrPeriod = new double[n][nPer];
    for( int q = 0; q < nPer; q++ ) {
      double[] e = solve( n, pa, pb, perShift[q], perQual[q], p.minQuality, refIdx, null );
      for( int i = 0; i < n; i++ ) r.nodeErrPeriod[i][q] = e[i];
    }
    r.nodeDrift = new double[n];
    for( int i = 0; i < n; i++ ) r.nodeDrift[i] = slopePerDay( r.periodTime, r.nodeErrPeriod[i] );
    r.nodeResid = new double[n];
    r.nodeSuspect = new boolean[n];
    trendResidual( nd, r.nodeErr, dt, r.nodeResid, r.nodeSuspect );
    return r;
  }

  //====================================================================
  private static double[] bandMask( int nfft, double df, double fmin, double fmax ) {
    double[] m = new double[nfft];
    double bw = Math.max( fmax - fmin, df );
    double edge = 0.2 * bw;
    for( int k = 0; k <= nfft / 2; k++ ) {
      double f = k * df, v;
      if( f < fmin - edge || f > fmax + edge ) v = 0.0;
      else if( f < fmin ) v = 0.5 - 0.5 * Math.cos( Math.PI * ( f - ( fmin - edge ) ) / edge );
      else if( f > fmax ) v = 0.5 - 0.5 * Math.cos( Math.PI * ( ( fmax + edge ) - f ) / edge );
      else v = 1.0;
      if( fmin - edge <= 0.0 && k == 0 ) v = 0.0;  // no DC
      m[k] = v;
      if( k > 0 && k < nfft / 2 ) m[nfft - k] = v;
    }
    return m;
  }
  private static double[] tukey( int n, double frac ) {
    double[] w = new double[n];
    int ne = Math.max( 1, (int)( frac * n ) );
    for( int i = 0; i < n; i++ ) {
      double v = 1.0;
      if( i < ne ) v = 0.5 - 0.5 * Math.cos( Math.PI * i / ne );
      else if( i >= n - ne ) v = 0.5 - 0.5 * Math.cos( Math.PI * ( n - 1 - i ) / ne );
      w[i] = v;
    }
    return w;
  }
  /** Window of nw samples starting at fractional index x0 of trace s -> pre-processed spectrum re/im */
  private static void preprocess( float[] s, double x0, int nw, double[] taper, double[] mask, Params p,
                                  double[] re, double[] im ) {
    int nfft = re.length;
    Arrays.fill( re, 0.0 );
    Arrays.fill( im, 0.0 );
    double mean = 0.0;
    for( int i = 0; i < nw; i++ ) {
      re[i] = interp( s, x0 + i );
      mean += re[i];
    }
    mean /= nw;
    for( int i = 0; i < nw; i++ ) re[i] = ( re[i] - mean ) * taper[i];
    fft( re, im );
    for( int f = 0; f < nfft; f++ ) { re[f] *= mask[f]; im[f] *= mask[f]; }
    if( p.oneBit ) {
      ifft( re, im );
      for( int i = 0; i < nfft; i++ ) {
        re[i] = ( i < nw ) ? Math.signum( re[i] ) * taper[i] : 0.0;
        im[i] = 0.0;
      }
      fft( re, im );
      for( int f = 0; f < nfft; f++ ) { re[f] *= mask[f]; im[f] *= mask[f]; }
    }
    if( p.whiten ) {
      // divide by the amplitude spectrum smoothed over ~ +/- 0.05 Hz (at least 3 bins)
      int h = nfft / 2 + 1;
      double[] amp = new double[h];
      for( int f = 0; f < h; f++ ) amp[f] = Math.hypot( re[f], im[f] );
      int half = Math.max( 1, nfft / 200 );
      double[] sm = movingAverage( amp, half );
      double eps = 0.0;
      for( double v : sm ) eps = Math.max( eps, v );
      eps *= 1.0e-6;
      for( int f = 0; f < h; f++ ) {
        double g = mask[f] / ( sm[f] + eps );
        re[f] *= g; im[f] *= g;
        if( f > 0 && f < nfft / 2 ) { re[nfft - f] *= g; im[nfft - f] *= g; }
      }
    }
  }
  private static double[] movingAverage( double[] a, int half ) {
    int n = a.length;
    double[] c = new double[n + 1];
    for( int i = 0; i < n; i++ ) c[i + 1] = c[i] + a[i];
    double[] o = new double[n];
    for( int i = 0; i < n; i++ ) {
      int i0 = Math.max( 0, i - half ), i1 = Math.min( n, i + half + 1 );
      o[i] = ( c[i1] - c[i0] ) / ( i1 - i0 );
    }
    return o;
  }

  /**
   * Symmetry measurement of a stacked correlation c[0..2nl] (index nl = lag 0).
   * @return { shift in SAMPLES (a = e_b - e_a), quality (normalized correlation of the two
   *           branches at the best shift), max|C| positive branch, max|C| negative branch }
   */
  public static double[] measure( double[] c, int nl, int kmin ) {
    int m = nl + 1;
    double[] P = new double[m], N = new double[m];
    double ap = 0.0, an = 0.0, ep = 0.0, en = 0.0;
    for( int k = 0; k < m; k++ ) {
      if( k < kmin ) continue;
      P[k] = c[nl + k];
      N[k] = c[nl - k];
      ap = Math.max( ap, Math.abs( P[k] ) );
      an = Math.max( an, Math.abs( N[k] ) );
      ep += P[k] * P[k];
      en += N[k] * N[k];
    }
    if( ep <= 0.0 || en <= 0.0 ) return new double[] { Double.NaN, 0.0, ap, an };
    // R(s) = sum_k P(k) N(k+s), via FFT
    int nfft = 1;
    while( nfft < 2 * m ) nfft <<= 1;
    double[] pr = new double[nfft], pi = new double[nfft], qr = new double[nfft], qi = new double[nfft];
    for( int k = 0; k < m; k++ ) { pr[k] = P[k]; qr[k] = N[k]; }
    fft( pr, pi );
    fft( qr, qi );
    for( int f = 0; f < nfft; f++ ) {
      double a = pr[f], b = pi[f], cr = qr[f], ci = qi[f];
      // conj(P) * Q
      pr[f] = a * cr + b * ci;
      pi[f] = a * ci - b * cr;
    }
    ifft( pr, pi );
    int smax = m - 1 - kmin;
    double best = Double.NEGATIVE_INFINITY;
    int sBest = 0;
    for( int sft = -smax; sft <= smax; sft++ ) {
      double v = pr[ sft >= 0 ? sft : nfft + sft ];
      if( v > best ) { best = v; sBest = sft; }
    }
    double frac = 0.0;
    if( sBest > -smax && sBest < smax ) {
      double y0 = pr[ ( sBest - 1 ) >= 0 ? sBest - 1 : nfft + sBest - 1 ];
      double y2 = pr[ ( sBest + 1 ) >= 0 ? sBest + 1 : nfft + sBest + 1 ];
      double den = y0 - 2.0 * best + y2;
      if( den < 0.0 ) frac = 0.5 * ( y0 - y2 ) / den;
    }
    double q = best / Math.sqrt( ep * en );
    // P(k) = N(k + 2a)  ->  best s = 2a
    return new double[] { 0.5 * ( sBest + frac ), q, ap, an };
  }

  /**
   * Least squares e_b - e_a = shift (weights = quality^2, pairs below minQuality skipped), with
   * e[ref] = 0 (ref >= 0) or sum(e) = 0. Nodes without any valid pair get NaN.
   */
  public static double[] solve( int n, int[] pa, int[] pb, double[] shift, double[] qual, double minQ,
                                int ref, int[] used ) {
    double[][] A = new double[n][n];
    double[] rhs = new double[n];
    int[] cnt = new int[n];
    for( int k = 0; k < pa.length; k++ ) {
      if( Double.isNaN( shift[k] ) || !( qual[k] >= minQ ) ) continue;
      int a = pa[k], b = pb[k];
      double w = qual[k] * qual[k];
      A[a][a] += w; A[b][b] += w; A[a][b] -= w; A[b][a] -= w;
      rhs[b] += w * shift[k];
      rhs[a] -= w * shift[k];
      cnt[a]++; cnt[b]++;
    }
    if( used != null ) System.arraycopy( cnt, 0, used, 0, n );
    double big = 0.0;
    for( int i = 0; i < n; i++ ) big = Math.max( big, A[i][i] );
    if( big <= 0.0 ) {
      double[] e = new double[n];
      Arrays.fill( e, Double.NaN );
      return e;
    }
    big *= 10.0;
    for( int i = 0; i < n; i++ ) {
      if( cnt[i] == 0 ) { A[i][i] += big; continue; }  // isolated node: pin to 0, reported NaN
      if( ref >= 0 ) {
        if( i == ref ) A[i][i] += big;
      }
      else {
        for( int j = 0; j < n; j++ ) if( cnt[j] > 0 ) A[i][j] += big / n;   // sum(e) = 0
      }
    }
    double[] e = gauss( A, rhs );
    for( int i = 0; i < n; i++ ) if( cnt[i] == 0 ) e[i] = Double.NaN;
    if( ref < 0 ) {   // remove any residual mean
      double m = 0.0; int c = 0;
      for( double v : e ) if( !Double.isNaN( v ) ) { m += v; c++; }
      if( c > 0 ) { m /= c; for( int i = 0; i < n; i++ ) if( !Double.isNaN( e[i] ) ) e[i] -= m; }
    }
    return e;
  }
  private static double[] gauss( double[][] A, double[] b ) {
    int n = b.length;
    double[][] M = new double[n][];
    for( int i = 0; i < n; i++ ) M[i] = A[i].clone();
    double[] x = b.clone();
    for( int c = 0; c < n; c++ ) {
      int piv = c;
      for( int r = c + 1; r < n; r++ ) if( Math.abs( M[r][c] ) > Math.abs( M[piv][c] ) ) piv = r;
      double[] t = M[c]; M[c] = M[piv]; M[piv] = t;
      double tb = x[c]; x[c] = x[piv]; x[piv] = tb;
      double d = M[c][c];
      if( Math.abs( d ) < 1.0e-300 ) continue;
      for( int r = c + 1; r < n; r++ ) {
        double f = M[r][c] / d;
        if( f == 0.0 ) continue;
        for( int k = c; k < n; k++ ) M[r][k] -= f * M[c][k];
        x[r] -= f * x[c];
      }
    }
    for( int r = n - 1; r >= 0; r-- ) {
      double sum = x[r];
      for( int k = r + 1; k < n; k++ ) sum -= M[r][k] * x[k];
      x[r] = ( Math.abs( M[r][r] ) < 1.0e-300 ) ? 0.0 : sum / M[r][r];
    }
    return x;
  }
  /**
   * Residual of the clock errors after removing a robust linear trend with the node number.
   * Independent clocks do not line up along the line: a smooth trend is usually a bias of the
   * symmetry method caused by directional (non-isotropic) noise, which grows with the distance;
   * a node with a real clock error sticks out of the trend. Fit: least squares, repeated without
   * the points beyond 3 robust standard deviations (MAD). Suspect: |residual| > max(3 sigma, 2 dt).
   */
  static void trendResidual( int[] x, double[] e, double dt, double[] resid, boolean[] suspect ) {
    int n = x.length;
    boolean[] use = new boolean[n];
    for( int i = 0; i < n; i++ ) use[i] = !Double.isNaN( e[i] );
    double a = 0.0, b = 0.0, sigma = 0.0;
    for( int iter = 0; iter < 4; iter++ ) {
      double sx = 0, sy = 0, sxx = 0, sxy = 0; int c = 0;
      for( int i = 0; i < n; i++ ) if( use[i] ) { sx += x[i]; sy += e[i]; sxx += (double)x[i] * x[i]; sxy += x[i] * e[i]; c++; }
      if( c < 3 ) { a = ( c > 0 ) ? sy / c : 0.0; b = 0.0; }
      else {
        double den = c * sxx - sx * sx;
        b = ( den != 0.0 ) ? ( c * sxy - sx * sy ) / den : 0.0;
        a = ( sy - b * sx ) / c;
      }
      double[] ar = new double[n]; int m = 0;
      for( int i = 0; i < n; i++ ) if( use[i] ) ar[m++] = Math.abs( e[i] - ( a + b * x[i] ) );
      double[] s = Arrays.copyOf( ar, m );
      Arrays.sort( s );
      sigma = ( m > 0 ) ? 1.4826 * s[m / 2] : 0.0;
      boolean changed = false;
      for( int i = 0; i < n; i++ ) {
        if( Double.isNaN( e[i] ) ) continue;
        boolean u = Math.abs( e[i] - ( a + b * x[i] ) ) <= Math.max( 3.0 * sigma, 2.0 * dt );
        if( u != use[i] ) { use[i] = u; changed = true; }
      }
      if( !changed ) break;
    }
    for( int i = 0; i < n; i++ ) {
      resid[i] = Double.isNaN( e[i] ) ? Double.NaN : e[i] - ( a + b * x[i] );
      suspect[i] = !Double.isNaN( resid[i] ) && Math.abs( resid[i] ) > Math.max( 3.0 * sigma, 2.0 * dt );
    }
  }

  /** Minimum time span between the first and last period for a drift to be reported [s] */
  public static final double MIN_DRIFT_SPAN_SEC = 3600.0;

  /**
   * Linear-regression slope of e (s) against time (s), in s/day; NaN with fewer than 2 points or
   * when the points span less than MIN_DRIFT_SPAN_SEC (too short to extrapolate to a day).
   */
  private static double slopePerDay( double[] t, double[] e ) {
    double st = 0, se = 0, stt = 0, ste = 0;
    double tmin = Double.POSITIVE_INFINITY, tmax = Double.NEGATIVE_INFINITY;
    int c = 0;
    for( int i = 0; i < t.length; i++ ) {
      if( Double.isNaN( e[i] ) ) continue;
      st += t[i]; se += e[i]; stt += t[i] * t[i]; ste += t[i] * e[i]; c++;
      tmin = Math.min( tmin, t[i] ); tmax = Math.max( tmax, t[i] );
    }
    if( c < 2 || tmax - tmin < MIN_DRIFT_SPAN_SEC ) return Double.NaN;
    double den = c * stt - st * st;
    if( den <= 0.0 ) return Double.NaN;
    return ( c * ste - st * se ) / den * 86400.0;
  }

  //====================================================================
  /** In-place radix-2 FFT (forward), length power of 2 */
  public static void fft( double[] re, double[] im ) {
    csPluginFX.fft( re, im );
  }
  /** In-place inverse FFT with 1/n scaling */
  public static void ifft( double[] re, double[] im ) {
    int n = re.length;
    for( int i = 0; i < n; i++ ) im[i] = -im[i];
    csPluginFX.fft( re, im );
    for( int i = 0; i < n; i++ ) { re[i] /= n; im[i] = -im[i] / n; }
  }
}
