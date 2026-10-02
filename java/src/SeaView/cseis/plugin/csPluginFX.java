/* SeaView plugin F-X spectrum: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import cseis.general.csUnits;
import cseis.jni.csVirtualSeismicReader;
import cseis.seis.csHeader;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicTraceBuffer;
import cseis.seis.csTraceBuffer;
import java.awt.GridLayout;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 * Plugin: F-X amplitude spectrum.<br>
 * Amplitude spectrum of every trace of the active pane, in a time window, opened in a new pane
 * (horizontal: traces, vertical: frequency [Hz]). Headers of the input traces are kept.
 */
public class csPluginFX implements csISeaViewPlugin {
  private static int ourCounter = 0;
  // Last used parameters (kept between runs)
  private String myTmin = "";
  private String myTmax = "";
  private String myFmax = "";
  private boolean myDB = true;
  private boolean myNormTrace = false;

  @Override
  public String getName() {
    return "F-X spectrum...";
  }
  @Override
  public String getDescription() {
    return "Amplitude spectrum of each trace (F-X) in a time window, shown in a new pane";
  }
  @Override
  public void run( csIPluginContext ctx ) {
    csISeismicTraceBuffer buffer = ctx.getTraceBuffer();
    if( buffer == null || buffer.numTraces() == 0 ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "No data in the active pane", "F-X spectrum", JOptionPane.WARNING_MESSAGE );
      return;
    }
    if( ctx.getVerticalDomain() == csUnits.DOMAIN_FREQ ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "Data in the active pane are already in the frequency domain", "F-X spectrum", JOptionPane.WARNING_MESSAGE );
      return;
    }
    float dt = ctx.getSampleInt();               // [ms]
    int ns = ctx.getNumSamples();
    float tend = (ns-1)*dt;
    String unit = ( ctx.getVerticalDomain() == csUnits.DOMAIN_DEPTH ) ? "m" : "ms";
    String funit = ( ctx.getVerticalDomain() == csUnits.DOMAIN_DEPTH ) ? "1/km" : "Hz";
    float fnyq = 500.0f / dt;

    //--- Parameter dialog
    JTextField textTmin = new JTextField( myTmin.length() > 0 ? myTmin : "0" );
    JTextField textTmax = new JTextField( myTmax.length() > 0 ? myTmax : ""+tend );
    JTextField textFmax = new JTextField( myFmax.length() > 0 ? myFmax : ""+fnyq );
    JCheckBox boxDB   = new JCheckBox( "Amplitude in dB", myDB );
    JCheckBox boxNorm = new JCheckBox( "Normalize each trace", myNormTrace );
    JPanel panel = new JPanel( new GridLayout(5,2,6,4) );
    panel.add( new JLabel("Window start [" + unit + "]:") ); panel.add( textTmin );
    panel.add( new JLabel("Window end [" + unit + "]:") );   panel.add( textTmax );
    panel.add( new JLabel("Max frequency [" + funit + "]:") ); panel.add( textFmax );
    panel.add( boxDB );   panel.add( new JLabel("") );
    panel.add( boxNorm ); panel.add( new JLabel("") );
    int option = JOptionPane.showConfirmDialog( ctx.getParentFrame(), panel, "F-X spectrum - " + ctx.getTitle(),
                                                JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE );
    if( option != JOptionPane.OK_OPTION ) return;
    float tmin, tmax, fmax;
    try {
      tmin = Float.parseFloat( textTmin.getText().trim() );
      tmax = Float.parseFloat( textTmax.getText().trim() );
      fmax = Float.parseFloat( textFmax.getText().trim() );
    }
    catch( NumberFormatException e ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "Invalid number", "F-X spectrum", JOptionPane.ERROR_MESSAGE );
      return;
    }
    myTmin = textTmin.getText().trim(); myTmax = textTmax.getText().trim(); myFmax = textFmax.getText().trim();
    myDB = boxDB.isSelected(); myNormTrace = boxNorm.isSelected();
    int i1 = Math.max( 0, Math.round( tmin / dt ) );
    int i2 = Math.min( ns-1, Math.round( tmax / dt ) );
    if( i2 - i1 + 1 < 8 ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "Time window too small (at least 8 samples)", "F-X spectrum", JOptionPane.ERROR_MESSAGE );
      return;
    }

    //--- Compute
    float[][] amp = computeFX( buffer, i1, i2, dt, myNormTrace, myDB );
    int nfAll = amp[0].length;
    float df = 1000.0f / ( dt * fftLength( i2-i1+1 ) );
    int nf = Math.max( 2, Math.min( nfAll, (int)Math.floor( fmax / df ) + 1 ) );

    csHeaderDef[] hdrDef = ctx.getHeaderDef();
    csVirtualSeismicReader reader = new csVirtualSeismicReader( nf, hdrDef.length, df, hdrDef, csUnits.DOMAIN_FREQ );
    csTraceBuffer out = reader.retrieveTraceBuffer();
    for( int itrc = 0; itrc < buffer.numTraces(); itrc++ ) {
      float[] s = new float[nf];
      System.arraycopy( amp[itrc], 0, s, 0, nf );
      csHeader[] hin = buffer.headerValues( itrc );
      csHeader[] hout = new csHeader[hdrDef.length];
      for( int ih = 0; ih < hdrDef.length; ih++ ) hout[ih] = ( ih < hin.length ) ? new csHeader( hin[ih] ) : new csHeader( 0 );
      out.addTrace( s, hout );
    }
    ourCounter++;
    String title = "FX" + ourCounter + " " + ctx.getTitle() + " (" + (int)(i1*dt) + "-" + (int)(i2*dt) + unit + ( myDB ? ", dB" : "" ) + ")";
    ctx.openNewPane( reader, title );
  }

  /** FFT length used for a window of n samples (power of 2, at least 2*n for a smooth spectrum) */
  public static int fftLength( int n ) {
    int nfft = 1;
    while( nfft < 2*n ) nfft *= 2;
    return nfft;
  }
  /**
   * Amplitude spectra of all traces in samples i1..i2 (Hann taper, zero padding).
   * @return amp[trace][frequency], frequency increment 1000/(dt*fftLength(i2-i1+1))
   */
  public static float[][] computeFX( csISeismicTraceBuffer buffer, int i1, int i2, float dt, boolean normTrace, boolean db ) {
    int n = i2 - i1 + 1;
    int nfft = fftLength( n );
    int nf = nfft/2 + 1;
    int ntr = buffer.numTraces();
    float[][] amp = new float[ntr][nf];
    double[] re = new double[nfft];
    double[] im = new double[nfft];
    double[] taper = new double[n];
    for( int i = 0; i < n; i++ ) taper[i] = 0.5 - 0.5*Math.cos( 2.0*Math.PI*(i+0.5)/n );
    double scale = dt * 0.001;     // amplitude in [U/Hz]
    float maxAll = 0.0f;
    for( int itrc = 0; itrc < ntr; itrc++ ) {
      float[] s = buffer.samples( itrc );
      java.util.Arrays.fill( re, 0.0 );
      java.util.Arrays.fill( im, 0.0 );
      double mean = 0.0;
      for( int i = 0; i < n; i++ ) mean += s[i1+i];
      mean /= n;
      for( int i = 0; i < n; i++ ) re[i] = ( s[i1+i] - mean ) * taper[i];
      fft( re, im );
      float maxTrace = 0.0f;
      for( int k = 0; k < nf; k++ ) {
        amp[itrc][k] = (float)( scale * Math.sqrt( re[k]*re[k] + im[k]*im[k] ) );
        maxTrace = Math.max( maxTrace, amp[itrc][k] );
      }
      if( normTrace && maxTrace > 0.0f ) {
        for( int k = 0; k < nf; k++ ) amp[itrc][k] /= maxTrace;
      }
      for( int k = 0; k < nf; k++ ) maxAll = Math.max( maxAll, amp[itrc][k] );
    }
    if( db ) {
      // dB relative to the maximum, clipped at -100 dB
      double ref = ( maxAll > 0.0f ) ? maxAll : 1.0;
      for( int itrc = 0; itrc < ntr; itrc++ ) {
        for( int k = 0; k < nf; k++ ) {
          double r = amp[itrc][k] / ref;
          amp[itrc][k] = (float)Math.max( -100.0, 20.0*Math.log10( Math.max( r, 1e-30 ) ) );
        }
      }
    }
    return amp;
  }
  /** In-place radix-2 complex FFT (forward, no scaling). Length must be a power of 2. */
  public static void fft( double[] re, double[] im ) {
    int n = re.length;
    for( int i = 1, j = 0; i < n; i++ ) {
      int bit = n >> 1;
      for( ; (j & bit) != 0; bit >>= 1 ) j ^= bit;
      j ^= bit;
      if( i < j ) {
        double t = re[i]; re[i] = re[j]; re[j] = t;
        t = im[i]; im[i] = im[j]; im[j] = t;
      }
    }
    for( int len = 2; len <= n; len <<= 1 ) {
      double ang = -2.0 * Math.PI / len;
      double wr = Math.cos( ang ), wi = Math.sin( ang );
      for( int i = 0; i < n; i += len ) {
        double cr = 1.0, ci = 0.0;
        for( int k = 0; k < len/2; k++ ) {
          int a = i + k, b = i + k + len/2;
          double xr = re[b]*cr - im[b]*ci;
          double xi = re[b]*ci + im[b]*cr;
          re[b] = re[a] - xr; im[b] = im[a] - xi;
          re[a] += xr;        im[a] += xi;
          double t = cr*wr - ci*wi;
          ci = cr*wi + ci*wr;
          cr = t;
        }
      }
    }
  }
}
