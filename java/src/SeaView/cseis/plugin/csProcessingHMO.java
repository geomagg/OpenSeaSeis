/* SeaView processing step HMO: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import cseis.processing.csIProcessing;
import cseis.seis.csDataBuffer;
import cseis.seis.csHeader;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicTraceBuffer;
import cseis.seis.csSeismicData;
import java.awt.GridLayout;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 * On-screen processing step: hyperbolic moveout correction with a chosen velocity.
 * <ul>
 * <li> HMO (shift): each trace is shifted by dt = (sqrt(x^2+z^2) - z)/v, without stretch.
 *      Flattens the direct arrival (and its mirror/multiples) of OBN receiver gathers, z = node depth.
 *      Optionally the event is flattened at a given time instead of at its zero-offset time.
 * <li> NMO: t0 = sqrt(t^2 - x^2/v^2), with optional stretch mute.
 * </ul>
 * x comes from an offset header, or from the source/receiver coordinates (sou_x/sou_y, rec_x/rec_y).
 * z is a constant or a trace header (absolute value used).
 */
public class csProcessingHMO implements csIProcessing {
  public static final String NAME = "HMO";
  private static final String[] MODES = { "HMO (shift, no stretch)", "NMO (hyperbolic stretch)" };

  private csHeaderDef[] myHdrDef;
  private float mySampleInt;
  private boolean myIsActive;
  private JPanel myPanel;
  private JComboBox<String> myComboMode;
  private JTextField myTextVel;
  private JTextField myTextOffset;
  private JTextField myTextDepth;
  private JTextField myTextFlat;
  private JTextField myTextMute;
  private JCheckBox myBoxInverse;

  // Parameters
  private int myMode;          // 0 HMO shift, 1 NMO
  private float myVel;         // [m/s]
  private int myHdrOffset;     // -1: from coordinates
  private int myHdrSx, myHdrSy, myHdrRx, myHdrRy;
  private float myDepth;
  private int myHdrDepth;      // -1: constant depth
  private float myFlatTime;    // [ms], <0: keep zero-offset time
  private float myMute;        // stretch mute [%], 0 = off
  private boolean myInverse;

  public csProcessingHMO( csHeaderDef[] hdrDef, float sampleInt ) {
    myHdrDef = hdrDef;
    mySampleInt = sampleInt;
    myIsActive = true;
    myComboMode = new JComboBox<String>( MODES );
    myTextVel    = new JTextField( "1500" );
    String offDefault = ( index("offset") >= 0 ) ? "offset" : "auto";
    myTextOffset = new JTextField( offDefault );
    myTextDepth  = new JTextField( "0" );
    myTextFlat   = new JTextField( "" );
    myTextMute   = new JTextField( "0" );
    myBoxInverse = new JCheckBox( "Inverse (undo correction)" );
    myTextVel.setToolTipText( "Velocity [m/s] (e.g. water velocity 1500)" );
    myTextOffset.setToolTipText( "<html>Offset header name, or <i>auto</i>: distance between (sou_x,sou_y) and (rec_x,rec_y)</html>" );
    myTextDepth.setToolTipText( "<html>HMO: depth of the node/source [m], a number or a header name (e.g. rec_z)<br>0 = linear-like moveout x/v</html>" );
    myTextFlat.setToolTipText( "<html>HMO: time [ms] where the corrected event is placed.<br>Empty = keep its zero-offset time z/v</html>" );
    myTextMute.setToolTipText( "NMO: stretch mute [%], 0 = no mute" );

    myPanel = new JPanel( new GridLayout(7,2,6,4) );
    myPanel.add( new JLabel("Mode:") );                     myPanel.add( myComboMode );
    myPanel.add( new JLabel("Velocity [m/s]:") );           myPanel.add( myTextVel );
    myPanel.add( new JLabel("Offset header (or auto):") );  myPanel.add( myTextOffset );
    myPanel.add( new JLabel("Depth z [m] or header:") );    myPanel.add( myTextDepth );
    myPanel.add( new JLabel("Flatten at time [ms]:") );     myPanel.add( myTextFlat );
    myPanel.add( new JLabel("NMO stretch mute [%]:") );     myPanel.add( myTextMute );
    myPanel.add( myBoxInverse );                            myPanel.add( new JLabel("") );
  }
  private int index( String name ) {
    if( myHdrDef == null ) return -1;
    for( int i = 0; i < myHdrDef.length; i++ ) {
      if( myHdrDef[i].name.compareTo( name ) == 0 ) return i;
    }
    return -1;
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
    myMode = myComboMode.getSelectedIndex();
    try {
      myVel = Float.parseFloat( myTextVel.getText().trim() );
    }
    catch( NumberFormatException e ) { return "Invalid velocity"; }
    if( !(myVel > 0.0f) ) return "Velocity must be > 0";
    String off = myTextOffset.getText().trim();
    myHdrSx = index("sou_x"); myHdrSy = index("sou_y"); myHdrRx = index("rec_x"); myHdrRy = index("rec_y");
    if( off.length() == 0 || off.equalsIgnoreCase("auto") ) {
      myHdrOffset = -1;
      if( myHdrSx < 0 || myHdrRx < 0 ) return "Offset 'auto' needs headers sou_x and rec_x. Give an offset header name instead.";
    }
    else {
      myHdrOffset = index( off );
      if( myHdrOffset < 0 ) return "Offset header '" + off + "' not found";
    }
    String dep = myTextDepth.getText().trim();
    myHdrDepth = -1;
    myDepth = 0.0f;
    if( dep.length() > 0 ) {
      try {
        myDepth = Math.abs( Float.parseFloat( dep ) );
      }
      catch( NumberFormatException e ) {
        myHdrDepth = index( dep );
        if( myHdrDepth < 0 ) return "Depth: '" + dep + "' is neither a number nor a header name";
      }
    }
    String flat = myTextFlat.getText().trim();
    myFlatTime = -1.0f;
    if( flat.length() > 0 ) {
      try { myFlatTime = Float.parseFloat( flat ); }
      catch( NumberFormatException e ) { return "Invalid flatten time"; }
    }
    try { myMute = Float.parseFloat( myTextMute.getText().trim() ); }
    catch( NumberFormatException e ) { return "Invalid stretch mute"; }
    myInverse = myBoxInverse.isSelected();
    return null;
  }

  private static double val( csHeader[] h, int i ) {
    return ( i >= 0 && i < h.length ) ? h[i].doubleValue() : 0.0;
  }
  @Override
  public void apply( csISeismicTraceBuffer in, csDataBuffer out ) {
    int ns = in.numSamples();
    float dt = mySampleInt;
    for( int itrc = 0; itrc < in.numTraces(); itrc++ ) {
      float[] s = in.samples( itrc );
      csHeader[] h = in.headerValues( itrc );
      double x;
      if( myHdrOffset >= 0 ) {
        x = Math.abs( val( h, myHdrOffset ) );
      }
      else {
        double dx = val( h, myHdrSx ) - val( h, myHdrRx );
        double dy = val( h, myHdrSy ) - val( h, myHdrRy );
        x = Math.sqrt( dx*dx + dy*dy );
      }
      float[] o = ( myMode == 0 ) ? hmoTrace( s, ns, dt, x, depthOf( h ) ) : nmoTrace( s, ns, dt, x );
      out.addDataTrace( new csSeismicData( o ) );
    }
  }
  private double depthOf( csHeader[] h ) {
    return ( myHdrDepth >= 0 ) ? Math.abs( val( h, myHdrDepth ) ) : myDepth;
  }
  /** Linear interpolation of trace s at fractional sample index */
  private static float interp( float[] s, int ns, double fidx ) {
    if( fidx < 0.0 || fidx > ns-1 ) return 0.0f;
    int i = (int)fidx;
    if( i >= ns-1 ) return s[ns-1];
    double w = fidx - i;
    return (float)( (1.0-w)*s[i] + w*s[i+1] );
  }
  /** HMO: time shift without stretch */
  public float[] hmoTrace( float[] s, int ns, float dt, double x, double z ) {
    double tEvent = 1000.0 * Math.sqrt( x*x + z*z ) / myVel;           // [ms]
    double tRef   = ( myFlatTime >= 0.0f ) ? myFlatTime : 1000.0 * z / myVel;
    double shift  = tEvent - tRef;                                      // [ms]
    if( myInverse ) shift = -shift;
    float[] o = new float[ns];
    for( int i = 0; i < ns; i++ ) o[i] = interp( s, ns, i + shift/dt );
    return o;
  }
  /** NMO: t0 = sqrt(t^2 - x^2/v^2), or inverse */
  public float[] nmoTrace( float[] s, int ns, float dt, double x ) {
    double tx = 1000.0 * x / myVel;   // [ms]
    float[] o = new float[ns];
    for( int i = 0; i < ns; i++ ) {
      double t = i * dt;
      double tin;
      if( !myInverse ) {
        tin = Math.sqrt( t*t + tx*tx );
        if( myMute > 0.0f && t > 0.0 && ( tin/t - 1.0 )*100.0 > myMute ) continue;
        if( t <= 0.0 && tx > 0.0 && myMute > 0.0f ) continue;
      }
      else {
        double a = t*t - tx*tx;
        if( a < 0.0 ) continue;
        tin = Math.sqrt( a );
      }
      o[i] = interp( s, ns, tin/dt );
    }
    return o;
  }
}
