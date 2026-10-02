/* SeaView plugin interface: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

/**
 * Tool acting on the active SeaView pane (F-X spectrum, HMO...).<br>
 * Tools are shown in the menu 'Plugins' by small adapter classes in cseis.plugins
 * (see cseis.plugins.csSpkFXPlugin), registered in cseis/resources/seaview_plugins.txt,
 * and run with SeaView.runPlugin(), which gives them a csIPluginContext.
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
