/* M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugins;

import cseis.plugin.csPluginAlignTime;
import cseis.plugin.csPluginClockXcorr;
import cseis.seaview.SeaView;
import cseis.seaview.plugin.csPluginContext;
import cseis.seaview.plugin.csSeaViewPlugin;

/**
 * Menu 'Plugins' entries for the OBN clock QC:
 * <ul>
 * <li>"Alinhar traços no tempo absoluto..." (cseis.plugin.csPluginAlignTime)</li>
 * <li>"Clock: correlação cruzada dos nodes..." (cseis.plugin.csPluginClockXcorr)</li>
 * </ul>
 * Registered in cseis/resources/seaview_plugins.txt.
 */
public class csSpkClockPlugin implements csSeaViewPlugin {
  @Override
  public String getName() {
    return "Clock dos nodes (alinhamento + correlação cruzada)";
  }
  @Override
  public void install( csPluginContext context ) {
    final SeaView seaview = context.getSeaView();
    final csPluginAlignTime align = new csPluginAlignTime();
    final csPluginClockXcorr xcorr = new csPluginClockXcorr();
    context.addMenuItem( align.getName(), () -> seaview.runPlugin( align ) ).setToolTipText( align.getDescription() );
    context.addMenuItem( xcorr.getName(), () -> seaview.runPlugin( xcorr ) ).setToolTipText( xcorr.getDescription() );
  }
}
