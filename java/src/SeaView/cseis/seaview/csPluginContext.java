/* SeaView plugin context implementation: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.seaview;

import cseis.plugin.csIPluginContext;
import cseis.processing.csIProcessing;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicReader;
import cseis.seis.csISeismicTraceBuffer;
import javax.swing.JFrame;

/**
 * Gives plugins access to the active seismic pane (bundle) and to SeaView.
 */
public class csPluginContext implements csIPluginContext {
  private SeaView mySeaView;
  private csSeisPaneBundle myBundle;

  public csPluginContext( SeaView seaview, csSeisPaneBundle bundle ) {
    mySeaView = seaview;
    myBundle  = bundle;
  }
  @Override
  public JFrame getParentFrame() { return mySeaView; }
  @Override
  public csISeismicTraceBuffer getTraceBuffer() { return myBundle.getSeismicTraceBuffer(); }
  @Override
  public csHeaderDef[] getHeaderDef() { return myBundle.getHeaderDef(); }
  @Override
  public int getHeaderIndex( String name ) {
    csHeaderDef[] def = myBundle.getHeaderDef();
    if( def == null ) return -1;
    for( int i = 0; i < def.length; i++ ) {
      if( def[i].name.compareTo( name ) == 0 ) return i;
    }
    return -1;
  }
  @Override
  public float getSampleInt() { return myBundle.getSampleInt(); }
  @Override
  public int getNumSamples() { return myBundle.getNumSamples(); }
  @Override
  public int getVerticalDomain() { return myBundle.getVerticalDomain(); }
  @Override
  public String getTitle() { return myBundle.getTitle(); }
  @Override
  public void openNewPane( csISeismicReader reader, String title ) {
    // SeaView reads one data set at a time and silently ignores a request made while a read is
    // still running (e.g. a plugin opening two panes in a row): wait until the previous read is done.
    if( !mySeaView.isReadProcessOngoing() ) {
      mySeaView.readData( reader, title, SeaView.FORMAT_CSEIS, true );
      return;
    }
    javax.swing.Timer timer = new javax.swing.Timer( 150, null );
    timer.addActionListener( e -> {
      if( mySeaView.isReadProcessOngoing() ) return;
      timer.stop();
      mySeaView.readData( reader, title, SeaView.FORMAT_CSEIS, true );
    } );
    timer.start();
  }
  @Override
  public void addProcessingStep( csIProcessing proc ) {
    myBundle.setProcessingStep( proc );
  }
}
