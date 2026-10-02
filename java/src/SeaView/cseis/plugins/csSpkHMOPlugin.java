/* M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugins;

import cseis.plugin.csPluginHMO;
import cseis.seaview.SeaView;
import cseis.seaview.plugin.csPluginContext;
import cseis.seaview.plugin.csSeaViewPlugin;

/**
 * Menu 'Plugins' entry: HMO / NMO (velocity) (hyperbolic moveout with a chosen velocity, applied on screen).
 * Registered in cseis/resources/seaview_plugins.txt. The computation is in cseis.plugin.csPluginHMO.
 */
public class csSpkHMOPlugin implements csSeaViewPlugin {
  @Override
  public String getName() {
    return "HMO / NMO (velocity)";
  }
  @Override
  public void install( csPluginContext context ) {
    final SeaView seaview = context.getSeaView();
    final csPluginHMO tool = new csPluginHMO();
    context.addMenuItem( tool.getName(), () -> seaview.runPlugin( tool ) )
           .setToolTipText( tool.getDescription() );
  }
}
