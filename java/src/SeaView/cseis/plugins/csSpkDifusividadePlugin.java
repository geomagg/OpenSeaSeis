/* M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugins;

import cseis.plugin.csPluginDifusividade;
import cseis.seaview.SeaView;
import cseis.seaview.plugin.csPluginContext;
import cseis.seaview.plugin.csSeaViewPlugin;

/**
 * Menu 'Plugins' entry: Difusividade do campo -- SPAC / Bessel J0, f-k, hourly symmetry and neighbour delays of the
 * raw noise in the active pane. Registered in cseis/resources/seaview_plugins.txt.
 * The computation is in cseis.plugin.csDifusividade.
 */
public class csSpkDifusividadePlugin implements csSeaViewPlugin {
  @Override
  public String getName() {
    return "Difusividade do campo";
  }
  @Override
  public void install( csPluginContext context ) {
    final SeaView seaview = context.getSeaView();
    final csPluginDifusividade tool = new csPluginDifusividade();
    context.addMenuItem( tool.getName(), () -> seaview.runPlugin( tool ) )
           .setToolTipText( tool.getDescription() );
  }
}
