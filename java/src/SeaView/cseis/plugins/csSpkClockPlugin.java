/* M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugins;

import cseis.plugin.csPluginClockXcorr;
import cseis.seaview.SeaView;
import cseis.seaview.plugin.csPluginContext;
import cseis.seaview.plugin.csSeaViewPlugin;

/**
 * Menu 'Plugins' entry for the OBN clock QC: "Clock: correlação cruzada dos nodes..." (cseis.plugin.csPluginClockXcorr).
 * "Alinhar traços no tempo absoluto..." is a separate entry (csSpkAlinharPlugin), so that other plugins can be listed
 * between them. Registered in cseis/resources/seaview_plugins.txt.
 */
public class csSpkClockPlugin implements csSeaViewPlugin {
  @Override
  public String getName() {
    return "Clock dos nodes (correlação cruzada)";
  }
  @Override
  public void install( csPluginContext context ) {
    final SeaView seaview = context.getSeaView();
    final csPluginClockXcorr xcorr = new csPluginClockXcorr();
    context.addMenuItem( xcorr.getName(), () -> seaview.runPlugin( xcorr ) ).setToolTipText( xcorr.getDescription() );
  }
}
