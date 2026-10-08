/* M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugins;

import cseis.plugin.csPluginInterferometria;
import cseis.seaview.SeaView;
import cseis.seaview.plugin.csPluginContext;
import cseis.seaview.plugin.csSeaViewPlugin;

/**
 * Menu 'Plugins' entry: Interferometria -- virtual shot gather (VSG) by cross-correlation with a reference trace.
 * Registered in cseis/resources/seaview_plugins.txt. The computation is in cseis.plugin.csPluginInterferometria.
 */
public class csSpkInterferometriaPlugin implements csSeaViewPlugin {
  @Override
  public String getName() {
    return "Interferometria";
  }
  @Override
  public void install( csPluginContext context ) {
    final SeaView seaview = context.getSeaView();
    final csPluginInterferometria tool = new csPluginInterferometria();
    context.addMenuItem( tool.getName(), () -> seaview.runPlugin( tool ) )
           .setToolTipText( tool.getDescription() );
  }
}
