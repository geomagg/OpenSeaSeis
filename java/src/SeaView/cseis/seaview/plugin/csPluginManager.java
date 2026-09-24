package cseis.seaview.plugin;

import cseis.seaview.SeaView;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import java.io.File;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.ServiceConfigurationError;
import java.util.ServiceLoader;

/**
 * Descobre e instala plugins.
 * <p>
 * Fontes, nesta ordem:
 * <ol>
 * <li>Plugins embutidos, listados em {@code cseis/resources/seaview_plugins.txt}.</li>
 * <li>Plugins no classpath registrados via {@code META-INF/services}.</li>
 * <li>Jars em {@code -Dseaview.plugins.dir=...} (várias pastas separadas por ':' ou ';').</li>
 * <li>Jars em {@code <pasta do aplicativo>/plugins}.</li>
 * <li>Jars em {@code ~/.seaview/plugins}.</li>
 * </ol>
 * Um plugin com erro é ignorado (e reportado) sem impedir os demais.
 */
public final class csPluginManager {
  private static final List<String> loadedNames = new ArrayList<>();
  private static final List<String> failures = new ArrayList<>();

  private csPluginManager() {}

  public static void installAll( SeaView seaview, JMenuBar menuBar ) {
    JMenu menu = new JMenu( "Plugins" );
    menu.setMnemonic( 'P' );
    menuBar.add( menu );

    ClassLoader loader = buildClassLoader();
    Map<String,csSeaViewPlugin> plugins = new LinkedHashMap<>();
    loadBuiltInPlugins( loader, plugins );
    try {
      for( csSeaViewPlugin p : ServiceLoader.load( csSeaViewPlugin.class, loader ) ) {
        plugins.putIfAbsent( p.getClass().getName(), p );
      }
    }
    catch( ServiceConfigurationError e ) {
      failures.add( e.getMessage() );
    }

    for( csSeaViewPlugin p : plugins.values() ) {
      String name = safeName( p );
      try {
        p.install( new csPluginContext( seaview, menu, name ) );
        loadedNames.add( name );
      }
      catch( Throwable t ) {
        t.printStackTrace();
        failures.add( name + ": " + t );
      }
    }

    if( menu.getItemCount() > 0 ) menu.addSeparator();
    JMenuItem about = new JMenuItem( "Plugins carregados..." );
    about.addActionListener( e -> JOptionPane.showMessageDialog( seaview, summary(), "Plugins", JOptionPane.INFORMATION_MESSAGE ) );
    menu.add( about );

    if( !failures.isEmpty() ) System.err.println( "SeaView: plugins com erro:\n  " + String.join( "\n  ", failures ) );
  }

  /** Plugins listados em cseis/resources/seaview_plugins.txt (compilados junto com o SeaView por make_java.sh). */
  private static void loadBuiltInPlugins( ClassLoader loader, Map<String,csSeaViewPlugin> plugins ) {
    java.io.InputStream in = csPluginManager.class.getResourceAsStream( "/cseis/resources/seaview_plugins.txt" );
    if( in == null ) return;
    try( java.io.BufferedReader r = new java.io.BufferedReader( new java.io.InputStreamReader( in, java.nio.charset.StandardCharsets.UTF_8 ) ) ) {
      String line;
      while( (line = r.readLine()) != null ) {
        line = line.trim();
        if( line.isEmpty() || line.startsWith( "#" ) ) continue;
        try {
          Object obj = Class.forName( line, true, loader ).getDeclaredConstructor().newInstance();
          plugins.putIfAbsent( line, (csSeaViewPlugin)obj );
        }
        catch( Throwable t ) {
          failures.add( line + ": " + t );
        }
      }
    }
    catch( java.io.IOException e ) {
      failures.add( "seaview_plugins.txt: " + e );
    }
  }

  public static List<String> getLoadedPluginNames() { return new ArrayList<>( loadedNames ); }

  private static String summary() {
    StringBuilder sb = new StringBuilder();
    sb.append( loadedNames.isEmpty() ? "Nenhum plugin carregado.\n" : "Carregados:\n" );
    for( String n : loadedNames ) sb.append( "  • " ).append( n ).append( '\n' );
    if( !failures.isEmpty() ) {
      sb.append( "\nCom erro:\n" );
      for( String f : failures ) sb.append( "  • " ).append( f ).append( '\n' );
    }
    sb.append( "\nPastas de plugins externos (.jar):\n" );
    for( File d : pluginDirs() ) sb.append( "  " ).append( d.getAbsolutePath() ).append( d.isDirectory() ? "" : "  (não existe)" ).append( '\n' );
    return sb.toString();
  }

  private static String safeName( csSeaViewPlugin p ) {
    try {
      String n = p.getName();
      return n == null || n.isEmpty() ? p.getClass().getSimpleName() : n;
    }
    catch( Throwable t ) {
      return p.getClass().getSimpleName();
    }
  }

  private static ClassLoader buildClassLoader() {
    ClassLoader parent = csPluginManager.class.getClassLoader();
    List<URL> jars = new ArrayList<>();
    for( File dir : pluginDirs() ) {
      File[] files = dir.listFiles( (d, n) -> n.toLowerCase().endsWith( ".jar" ) );
      if( files == null ) continue;
      for( File f : files ) {
        try {
          jars.add( f.toURI().toURL() );
        }
        catch( MalformedURLException e ) {
          failures.add( f + ": " + e );
        }
      }
    }
    return jars.isEmpty() ? parent : new URLClassLoader( jars.toArray( new URL[0] ), parent );
  }

  static List<File> pluginDirs() {
    List<File> dirs = new ArrayList<>();
    String prop = System.getProperty( "seaview.plugins.dir" );
    if( prop != null && !prop.isEmpty() ) {
      for( String s : prop.split( File.pathSeparator ) ) {
        if( !s.isEmpty() ) dirs.add( new File( s ) );
      }
    }
    File appDir = applicationDir();
    if( appDir != null ) dirs.add( new File( appDir, "plugins" ) );
    dirs.add( new File( System.getProperty( "user.home" ), ".seaview" + File.separator + "plugins" ) );
    return dirs;
  }

  /** Pasta raiz do aplicativo empacotado (jpackage) ou pasta do jar. */
  private static File applicationDir() {
    String appPath = System.getProperty( "jpackage.app-path" );  // .../SeaView/bin/SeaView
    if( appPath != null ) {
      File bin = new File( appPath ).getAbsoluteFile().getParentFile();
      return bin == null ? null : bin.getParentFile();
    }
    try {
      File jar = new File( csPluginManager.class.getProtectionDomain().getCodeSource().getLocation().toURI() );
      return jar.isFile() ? jar.getParentFile() : null;
    }
    catch( Exception e ) {
      return null;
    }
  }
}
