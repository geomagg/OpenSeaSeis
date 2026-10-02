/* M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugins;

import cseis.plugin.csPluginFX;
import cseis.seaview.SeaView;
import cseis.seaview.plugin.csPluginContext;
import cseis.seaview.plugin.csSeaViewPlugin;

/**
 * Menu 'Plugins' entry: F-X spectrum (amplitude spectrum of each trace in a time window, in a new pane).
 * Registered in cseis/resources/seaview_plugins.txt. The computation is in cseis.plugin.csPluginFX.
 */
public class csSpkFXPlugin implements csSeaViewPlugin {
  @Override
  public String getName() {
    return "F-X spectrum";
  }
  @Override
  public void install( csPluginContext context ) {
    final SeaView seaview = context.getSeaView();
    final csPluginFX tool = new csPluginFX();
    context.addMenuItem( tool.getName(), () -> seaview.runPlugin( tool ) )
           .setToolTipText( tool.getDescription() );
  }
}
