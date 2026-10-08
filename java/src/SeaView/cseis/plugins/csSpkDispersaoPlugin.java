/* M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugins;

import cseis.plugin.csPluginDispersao;
import cseis.seaview.SeaView;
import cseis.seaview.plugin.csPluginContext;
import cseis.seaview.plugin.csSeaViewPlugin;

/**
 * Menu 'Plugins' entry: Dispersao f-c -- phase velocity vs frequency image of the active pane (e.g. a VSG).
 * Registered in cseis/resources/seaview_plugins.txt. The computation is in cseis.plugin.csPluginDispersao.
 */
public class csSpkDispersaoPlugin implements csSeaViewPlugin {
  @Override
  public String getName() {
    return "Dispersao f-c";
  }
  @Override
  public void install( csPluginContext context ) {
    final SeaView seaview = context.getSeaView();
    final csPluginDispersao tool = new csPluginDispersao();
    context.addMenuItem( tool.getName(), () -> seaview.runPlugin( tool ) )
           .setToolTipText( tool.getDescription() );
  }
}
