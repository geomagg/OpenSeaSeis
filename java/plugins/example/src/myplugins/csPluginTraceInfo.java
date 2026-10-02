/* Example of an external SeaView plugin */

package myplugins;

import cseis.plugin.csIPluginContext;
import cseis.plugin.csISeaViewPlugin;
import cseis.seis.csISeismicTraceBuffer;
import javax.swing.JOptionPane;

/**
 * Example plugin: shows number of traces, samples, sample interval and amplitude range of the active pane.
 * Copy this directory to start a new plugin.
 */
public class csPluginTraceInfo implements csISeaViewPlugin {
  @Override
  public String getName() {
    return "Trace info (example plugin)";
  }
  @Override
  public String getDescription() {
    return "Example of an external plugin: basic information about the active pane";
  }
  @Override
  public void run( csIPluginContext ctx ) {
    csISeismicTraceBuffer buf = ctx.getTraceBuffer();
    float amin = Float.MAX_VALUE, amax = -Float.MAX_VALUE;
    for( int i = 0; i < buf.numTraces(); i++ ) {
      for( float s : buf.samples(i) ) { amin = Math.min( amin, s ); amax = Math.max( amax, s ); }
    }
    JOptionPane.showMessageDialog( ctx.getParentFrame(),
        ctx.getTitle() + "\n" +
        "Traces: " + buf.numTraces() + "\n" +
        "Samples: " + ctx.getNumSamples() + "  (interval " + ctx.getSampleInt() + ")\n" +
        "Trace headers: " + ctx.getHeaderDef().length + "\n" +
        "Amplitude: " + amin + " ... " + amax,
        getName(), JOptionPane.INFORMATION_MESSAGE );
  }
}
