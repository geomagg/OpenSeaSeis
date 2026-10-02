/* SeaView plugin manager: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import java.io.File;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.ServiceConfigurationError;
import java.util.ServiceLoader;

/**
 * Builds the list of SeaView plugins: the built-in ones (F-X spectrum, HMO) and the external ones
 * found in jar files in the plugin directories:
 * <ul>
 *  <li> 'plugins' next to SeaView.jar (e.g. lib_v3.00/lib/plugins)
 *  <li> ~/.seaview/plugins
 *  <li> directory given in the environment variable SEAVIEW_PLUGINS
 * </ul>
 */
public class csPluginManager {
  private List<csISeaViewPlugin> myPlugins;
  private List<String> myMessages;

  public csPluginManager() {
    myPlugins  = new ArrayList<csISeaViewPlugin>();
    myMessages = new ArrayList<String>();
    // Built-in plugins
    myPlugins.add( new csPluginFX() );
    myPlugins.add( new csPluginHMO() );
    loadExternal();
  }
  public List<csISeaViewPlugin> getPlugins() {
    return myPlugins;
  }
  /** @return Messages from loading external plugins (errors, loaded jars) */
  public List<String> getMessages() {
    return myMessages;
  }
  /** @return Plugin directories that are searched */
  public static List<File> pluginDirectories() {
    List<File> dirs = new ArrayList<File>();
    try {
      URL loc = csPluginManager.class.getProtectionDomain().getCodeSource().getLocation();
      File jarFile = new File( loc.toURI() );
      dirs.add( new File( jarFile.getParentFile(), "plugins" ) );
    }
    catch( Exception e ) {
      // ignore: location unknown
    }
    dirs.add( new File( System.getProperty("user.home"), ".seaview" + File.separator + "plugins" ) );
    String env = System.getenv( "SEAVIEW_PLUGINS" );
    if( env != null && env.length() > 0 ) dirs.add( new File(env) );
    return dirs;
  }
  private void loadExternal() {
    List<URL> urls = new ArrayList<URL>();
    for( File dir : pluginDirectories() ) {
      File[] files = dir.listFiles();
      if( files == null ) continue;
      for( File f : files ) {
        if( f.isFile() && f.getName().toLowerCase().endsWith(".jar") ) {
          try {
            urls.add( f.toURI().toURL() );
            myMessages.add( "Plugin jar: " + f.getPath() );
          }
          catch( Exception e ) {
            myMessages.add( "Cannot use plugin jar " + f.getPath() + ": " + e.getMessage() );
          }
        }
      }
    }
    if( urls.isEmpty() ) return;
    URLClassLoader loader = new URLClassLoader( urls.toArray( new URL[0] ), csPluginManager.class.getClassLoader() );
    ServiceLoader<csISeaViewPlugin> service = ServiceLoader.load( csISeaViewPlugin.class, loader );
    Iterator<csISeaViewPlugin> it = service.iterator();
    while( true ) {
      try {
        if( !it.hasNext() ) break;
        csISeaViewPlugin plugin = it.next();
        myPlugins.add( plugin );
        myMessages.add( "Plugin loaded: " + plugin.getName() + " (" + plugin.getClass().getName() + ")" );
      }
      catch( ServiceConfigurationError e ) {
        myMessages.add( "Error loading plugin: " + e.getMessage() );
      }
      catch( Throwable e ) {
        myMessages.add( "Error loading plugin: " + e );
        break;
      }
    }
  }
}
