/* SeaView plugin HMO: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import javax.swing.JOptionPane;

/**
 * Plugin: HMO / NMO with a chosen velocity, applied on screen to the active pane.<br>
 * Opens the processing dialog: change velocity (or depth, mode) and press Apply as often as needed.
 * Remove it with 'Active pane / Processing / Clear processing'.
 */
public class csPluginHMO implements csISeaViewPlugin {
  @Override
  public String getName() {
    return "HMO / NMO (velocity)...";
  }
  @Override
  public String getDescription() {
    return "Hyperbolic moveout correction with a chosen velocity (HMO shift for OBN direct arrival, or NMO), on screen";
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
