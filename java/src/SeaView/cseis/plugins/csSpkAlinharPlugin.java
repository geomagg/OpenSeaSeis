/* M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugins;

import cseis.plugin.csPluginAlignTime;
import cseis.seaview.SeaView;
import cseis.seaview.plugin.csPluginContext;
import cseis.seaview.plugin.csSeaViewPlugin;

/**
 * Menu 'Plugins' entry "Alinhar traços no tempo absoluto..." (cseis.plugin.csPluginAlignTime).
 * Registered in cseis/resources/seaview_plugins.txt.
 */
public class csSpkAlinharPlugin implements csSeaViewPlugin {
  @Override
  public String getName() {
    return "Alinhar traços no tempo absoluto";
  }
  @Override
  public void install( csPluginContext context ) {
    final SeaView seaview = context.getSeaView();
    final csPluginAlignTime align = new csPluginAlignTime();
    context.addMenuItem( align.getName(), () -> seaview.runPlugin( align ) ).setToolTipText( align.getDescription() );
  }
}
