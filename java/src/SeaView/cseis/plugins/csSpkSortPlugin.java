/* M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugins;

import cseis.plugin.csPluginSort;
import cseis.seaview.SeaView;
import cseis.seaview.plugin.csPluginContext;
import cseis.seaview.plugin.csSeaViewPlugin;

/**
 * Menu 'Plugins' entry: sort the traces of the active pane by one or two header keys (new pane).
 * Registered in cseis/resources/seaview_plugins.txt. The sort is in cseis.plugin.csPluginSort.
 */
public class csSpkSortPlugin implements csSeaViewPlugin {
  @Override
  public String getName() {
    return "Sort by header";
  }
  @Override
  public void install( csPluginContext context ) {
    final SeaView seaview = context.getSeaView();
    final csPluginSort tool = new csPluginSort();
    context.addMenuItem( tool.getName(), () -> seaview.runPlugin( tool ) )
           .setToolTipText( tool.getDescription() );
  }
}
