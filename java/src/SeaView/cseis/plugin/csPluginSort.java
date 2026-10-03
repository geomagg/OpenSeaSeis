/* SeaView plugin SORT: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import cseis.jni.csJNIDef;
import cseis.jni.csVirtualSeismicReader;
import cseis.seis.csHeader;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicTraceBuffer;
import cseis.seis.csTraceBuffer;
import java.awt.GridLayout;
import java.util.Arrays;
import java.util.Comparator;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

/**
 * Plugin: sorts the traces of the active pane by one or two trace header keys and shows the
 * result in a new pane (samples and headers move together; the original pane is not changed).
 * The sort is stable: traces with equal keys keep their original order.
 */
public class csPluginSort implements csISeaViewPlugin {
  private static final String NONE = "(none)";
  private static final String[] ORDERS = { "increasing", "decreasing" };
  private static int ourCounter = 0;
  // Last choices (kept between calls)
  private String myKey1 = null, myKey2 = NONE;
  private int myOrder1 = 0, myOrder2 = 0;
  private boolean myAbs1 = false, myAbs2 = false;

  @Override
  public String getName() {
    return "Sort by header...";
  }
  @Override
  public String getDescription() {
    return "Sort the traces of the active pane by one or two trace header keys, shown in a new pane";
  }

  @Override
  public void run( csIPluginContext ctx ) {
    csISeismicTraceBuffer buffer = ctx.getTraceBuffer();
    csHeaderDef[] hdrDef = ctx.getHeaderDef();
    if( buffer == null || buffer.numTraces() == 0 || hdrDef == null || hdrDef.length == 0 ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "No data in the active pane", "Sort", JOptionPane.WARNING_MESSAGE );
      return;
    }
    String[] names = new String[hdrDef.length];
    for( int i = 0; i < hdrDef.length; i++ ) names[i] = hdrDef[i].name;
    String[] names2 = new String[hdrDef.length+1];
    names2[0] = NONE;
    System.arraycopy( names, 0, names2, 1, names.length );

    JComboBox<String> comboKey1 = new JComboBox<String>( names );
    JComboBox<String> comboKey2 = new JComboBox<String>( names2 );
    JComboBox<String> comboOrd1 = new JComboBox<String>( ORDERS );
    JComboBox<String> comboOrd2 = new JComboBox<String>( ORDERS );
    JCheckBox boxAbs1 = new JCheckBox( "absolute value", myAbs1 );
    JCheckBox boxAbs2 = new JCheckBox( "absolute value", myAbs2 );
    String key1 = myKey1;
    if( key1 == null || ctx.getHeaderIndex( key1 ) < 0 ) key1 = ( ctx.getHeaderIndex("offset") >= 0 ) ? "offset" : names[0];
    comboKey1.setSelectedItem( key1 );
    comboKey2.setSelectedItem( ctx.getHeaderIndex( myKey2 ) >= 0 ? myKey2 : NONE );
    comboOrd1.setSelectedIndex( myOrder1 );
    comboOrd2.setSelectedIndex( myOrder2 );
    comboKey1.setMaximumRowCount( 20 );
    comboKey2.setMaximumRowCount( 20 );

    JPanel panel = new JPanel( new GridLayout(2,4,6,4) );
    panel.add( new JLabel("Primary key:") );   panel.add( comboKey1 ); panel.add( comboOrd1 ); panel.add( boxAbs1 );
    panel.add( new JLabel("Secondary key:") ); panel.add( comboKey2 ); panel.add( comboOrd2 ); panel.add( boxAbs2 );
    int option = JOptionPane.showConfirmDialog( ctx.getParentFrame(), panel, "Sort - " + ctx.getTitle(),
                                                JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE );
    if( option != JOptionPane.OK_OPTION ) return;

    myKey1 = (String)comboKey1.getSelectedItem();
    myKey2 = (String)comboKey2.getSelectedItem();
    myOrder1 = comboOrd1.getSelectedIndex();
    myOrder2 = comboOrd2.getSelectedIndex();
    myAbs1 = boxAbs1.isSelected();
    myAbs2 = boxAbs2.isSelected();
    int idx1 = ctx.getHeaderIndex( myKey1 );
    int idx2 = NONE.equals( myKey2 ) ? -1 : ctx.getHeaderIndex( myKey2 );

    int[] order = sortOrder( buffer, idx1, myOrder1 == 1, myAbs1, idx2, myOrder2 == 1, myAbs2 );

    csVirtualSeismicReader reader = new csVirtualSeismicReader( buffer.numSamples(), hdrDef.length, ctx.getSampleInt(), hdrDef, ctx.getVerticalDomain() );
    csTraceBuffer out = reader.retrieveTraceBuffer();
    for( int k = 0; k < order.length; k++ ) {
      int itrc = order[k];
      float[] s = buffer.samples( itrc ).clone();
      csHeader[] hin = buffer.headerValues( itrc );
      csHeader[] hout = new csHeader[hdrDef.length];
      for( int ih = 0; ih < hdrDef.length; ih++ ) hout[ih] = ( ih < hin.length ) ? new csHeader( hin[ih] ) : new csHeader( 0 );
      out.addTrace( s, hout );
    }
    ourCounter++;
    String title = "SORT" + ourCounter + " " + ctx.getTitle() + " (" + describe( myKey1, myOrder1, myAbs1 )
                   + ( idx2 >= 0 ? ", " + describe( myKey2, myOrder2, myAbs2 ) : "" ) + ")";
    ctx.openNewPane( reader, title );
  }

  private static String describe( String key, int order, boolean abs ) {
    return ( abs ? "|" + key + "|" : key ) + ( order == 1 ? " desc" : "" );
  }

  /**
   * Trace order after sorting (stable).
   * @param idx1 header index of the primary key
   * @param idx2 header index of the secondary key, or -1
   * @return original trace indices in the new order
   */
  public static int[] sortOrder( csISeismicTraceBuffer buffer, final int idx1, final boolean desc1, final boolean abs1,
                                 final int idx2, final boolean desc2, final boolean abs2 ) {
    int ntr = buffer.numTraces();
    final csHeader[][] h = new csHeader[ntr][];
    Integer[] idx = new Integer[ntr];
    for( int i = 0; i < ntr; i++ ) {
      h[i] = buffer.headerValues( i );
      idx[i] = i;
    }
    Comparator<Integer> cmp = new Comparator<Integer>() {
      @Override
      public int compare( Integer a, Integer b ) {
        int c = compareKey( h[a], h[b], idx1, abs1 );
        if( desc1 ) c = -c;
        if( c == 0 && idx2 >= 0 ) {
          c = compareKey( h[a], h[b], idx2, abs2 );
          if( desc2 ) c = -c;
        }
        return c;
      }
    };
    Arrays.sort( idx, cmp );   // merge sort: stable
    int[] order = new int[ntr];
    for( int i = 0; i < ntr; i++ ) order[i] = idx[i];
    return order;
  }

  private static int compareKey( csHeader[] ha, csHeader[] hb, int ih, boolean abs ) {
    csHeader a = ( ih < ha.length ) ? ha[ih] : null;
    csHeader b = ( ih < hb.length ) ? hb[ih] : null;
    if( a == null || b == null ) return ( a == null ? ( b == null ? 0 : -1 ) : 1 );
    if( a.type() == csJNIDef.TYPE_STRING || b.type() == csJNIDef.TYPE_STRING ) {
      return String.valueOf( a.stringValue() ).compareTo( String.valueOf( b.stringValue() ) );
    }
    double va = a.doubleValue(), vb = b.doubleValue();
    if( abs ) { va = Math.abs( va ); vb = Math.abs( vb ); }
    return Double.compare( va, vb );
  }
}
