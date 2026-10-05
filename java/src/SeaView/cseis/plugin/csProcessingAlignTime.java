/* SeaView processing step: absolute-time alignment. M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import cseis.processing.csIProcessing;
import cseis.seis.csDataBuffer;
import cseis.seis.csHeader;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicTraceBuffer;
import cseis.seis.csSeismicData;
import java.awt.BorderLayout;
import java.time.Instant;
import java.util.Locale;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;

/**
 * On-screen processing step: aligns all traces in ABSOLUTE time, so that sample 0 of every
 * trace is the same instant -- the LATEST start time among the traces (trace headers time_samp1,
 * or time_year/day/hour/min/sec, see csClockXcorr.startTimeSec). Optionally every trace is also
 * cut where the EARLIEST-ending trace ends (samples after that are zeroed).<br>
 * Example: a trace starting at 7 s and one starting at 12 s both start at 12 s, and both end
 * where the 7 s trace ends.<br>
 * Remove with 'Active pane / Processing / Clear processing'.
 */
public class csProcessingAlignTime implements csIProcessing {
  public static final String NAME = "Align time";
  private final csHeaderDef[] myHdrDef;
  private final float mySampleInt;   // [ms]
  private boolean myIsActive = true;
  private final JPanel myPanel;
  private final JCheckBox myBoxCut;
  private boolean myCut = true;

  public csProcessingAlignTime( csHeaderDef[] hdrDef, float sampleInt, csISeismicTraceBuffer buffer ) {
    myHdrDef = hdrDef;
    mySampleInt = sampleInt;
    myBoxCut = new JCheckBox( "Cortar no fim comum (zera depois do fim do traço que termina primeiro)", true );
    myPanel = new JPanel( new BorderLayout( 4, 6 ) );
    myPanel.add( new JLabel( summary( buffer ) ), BorderLayout.CENTER );
    myPanel.add( myBoxCut, BorderLayout.SOUTH );
  }

  private String summary( csISeismicTraceBuffer buf ) {
    if( buf == null || buf.numTraces() == 0 ) return "Sem traços.";
    int n = buf.numTraces();
    double[] t0 = new double[n];
    int nMissing = 0;
    for( int i = 0; i < n; i++ ) {
      t0[i] = csClockXcorr.startTimeSec( buf.headerValues( i ), myHdrDef );
      if( Double.isNaN( t0[i] ) ) nMissing++;
    }
    if( nMissing == n ) return "<html>Nenhum traço tem tempo no header (time_samp1 ou time_day/hour/min/sec).</html>";
    double dt = mySampleInt / 1000.0;
    double tS = Double.NEGATIVE_INFINITY, tE = Double.POSITIVE_INFINITY, tMin = Double.POSITIVE_INFINITY;
    for( double t : t0 ) {
      if( Double.isNaN( t ) ) continue;
      tS = Math.max( tS, t );
      tMin = Math.min( tMin, t );
      tE = Math.min( tE, t + ( buf.numSamples() - 1 ) * dt );
    }
    StringBuilder sb = new StringBuilder( "<html>" );
    sb.append( String.format( Locale.US, "Início comum (o mais tardio): <b>%s</b><br>", fmt( tS ) ) );
    sb.append( String.format( Locale.US, "Fim comum (o mais cedo): <b>%s</b><br>", fmt( tE ) ) );
    sb.append( String.format( Locale.US, "Janela comum: %.3f s &nbsp; (maior deslocamento: %.3f s)<br>", tE - tS, tS - tMin ) );
    if( tE <= tS ) sb.append( "<font color=red>Os traços não se sobrepõem no tempo!</font><br>" );
    if( nMissing > 0 ) sb.append( String.format( Locale.US, "<font color=red>%d traço(s) sem tempo no header ficam como estão.</font><br>", nMissing ) );
    sb.append( "O tempo 0 da tela passa a ser o início comum.</html>" );
    return sb.toString();
  }
  static String fmt( double t ) {
    if( Double.isNaN( t ) || Double.isInfinite( t ) ) return "-";
    long ms = Math.round( t * 1000.0 );
    return Instant.ofEpochMilli( ms ).toString().replace( 'T', ' ' ).replace( "Z", " UTC" );
  }

  @Override
  public boolean isActive() { return myIsActive; }
  @Override
  public void setActive( boolean doSet ) { myIsActive = doSet; }
  @Override
  public String getName() { return NAME; }
  @Override
  public JPanel getParameterPanel() { return myPanel; }
  @Override
  public String retrieveParameters() {
    myCut = myBoxCut.isSelected();
    return null;
  }

  @Override
  public void apply( csISeismicTraceBuffer in, csDataBuffer out ) {
    int n = in.numTraces();
    int ns = in.numSamples();
    double dt = mySampleInt / 1000.0;
    double[] t0 = new double[n];
    double tS = Double.NEGATIVE_INFINITY, tE = Double.POSITIVE_INFINITY;
    for( int i = 0; i < n; i++ ) {
      csHeader[] h = in.headerValues( i );
      t0[i] = csClockXcorr.startTimeSec( h, myHdrDef );
      if( Double.isNaN( t0[i] ) ) continue;
      tS = Math.max( tS, t0[i] );
      tE = Math.min( tE, t0[i] + ( ns - 1 ) * dt );
    }
    int nKeep = ns;
    if( myCut && tE > tS ) nKeep = Math.min( ns, (int)Math.floor( ( tE - tS ) / dt + 1.0e-6 ) + 1 );
    for( int i = 0; i < n; i++ ) {
      float[] s = in.samples( i );
      float[] o = new float[ns];
      if( Double.isNaN( t0[i] ) || Double.isInfinite( tS ) ) {
        System.arraycopy( s, 0, o, 0, ns );
      }
      else {
        double shift = ( tS - t0[i] ) / dt;
        for( int k = 0; k < nKeep; k++ ) o[k] = (float)csClockXcorr.interp( s, k + shift );
      }
      out.addDataTrace( new csSeismicData( o ) );
    }
  }
}
