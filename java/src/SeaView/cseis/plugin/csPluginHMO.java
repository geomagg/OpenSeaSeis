/* SeaView plugin HMO: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import javax.swing.JOptionPane;

/**
 * Plugin: LMO / HMO with a chosen velocity, applied on screen to the active pane.<br>
 * Opens the processing dialog: change velocity (or depth, mode) and press Apply as often as needed.
 * Remove it with 'Active pane / Processing / Clear processing'.
 */
public class csPluginHMO implements csISeaViewPlugin {
  @Override
  public String getName() {
    return "LMO / HMO (velocity)...";
  }
  @Override
  public String getDescription() {
    return "Moveout with a chosen velocity, on screen: HMO (hyperbolic, t0 = sqrt(t^2-x^2/v^2)) or LMO (shift of each trace, flattens the direct arrival)";
  }
  @Override
  public void run( csIPluginContext ctx ) {
    if( ctx.getTraceBuffer() == null || ctx.getTraceBuffer().numTraces() == 0 ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "No data in the active pane", "HMO", JOptionPane.WARNING_MESSAGE );
      return;
    }
    ctx.addProcessingStep( new csProcessingHMO( ctx.getHeaderDef(), ctx.getSampleInt() ) );
  }
}
