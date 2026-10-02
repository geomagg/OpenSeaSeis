/* SeaView plugin interface: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

/**
 * SeaView plugin.<br>
 * A plugin appears as one item in the SeaView menu 'Plugins' and acts on the active data pane.<br>
 * External plugins: put a jar file in the 'plugins' directory next to SeaView.jar (or in ~/.seaview/plugins).
 * The jar must contain the plugin class(es) and the file META-INF/services/cseis.plugin.csISeaViewPlugin
 * listing the full class name of each plugin, one per line. The class needs a public no-argument constructor.
 */
public interface csISeaViewPlugin {
  /** @return Text of the menu item */
  public String getName();
  /** @return Short description (menu tooltip) */
  public String getDescription();
  /**
   * Run the plugin on the active data pane
   * @param context Access to the active data pane and to SeaView
   */
  public void run( csIPluginContext context );
}
