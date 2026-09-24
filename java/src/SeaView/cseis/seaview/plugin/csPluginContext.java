package cseis.seaview.plugin;

import cseis.seaview.SeaView;
import cseis.seaview.csSeisPaneBundle;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicTraceBuffer;

import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.KeyStroke;
import java.util.function.Consumer;

/**
 * O que um plugin pode acessar no SeaView: janela, painel ativo, dados e menu.
 */
public final class csPluginContext {
  private final SeaView mySeaView;
  private final JMenu myMenu;
  private final String myPluginName;

  csPluginContext( SeaView seaview, JMenu menu, String pluginName ) {
    mySeaView = seaview;
    myMenu = menu;
    myPluginName = pluginName;
  }

  /** Janela principal (use como "parent" de diálogos). */
  public SeaView getSeaView() { return mySeaView; }

  /** Menu "Plugins" (para submenus ou itens personalizados). */
  public JMenu getPluginMenu() { return myMenu; }

  /** Painel sísmico ativo, ou null se não houver arquivo aberto. */
  public csSeisPaneBundle getActiveBundle() { return mySeaView.getActiveBundle(); }

  /** Traços do painel ativo, ou null. */
  public csISeismicTraceBuffer getActiveTraces() {
    csSeisPaneBundle b = getActiveBundle();
    return b == null ? null : b.getTraceBuffer();
  }

  /** Cabeçalhos do painel ativo (vazio se não houver painel). */
  public csHeaderDef[] getActiveHeaderDefs() {
    csSeisPaneBundle b = getActiveBundle();
    return b == null ? new csHeaderDef[0] : b.getTraceHeaderDefs();
  }

  /** Abre um arquivo sísmico (SEG-Y, CSeis, SU, SEG-D, RSF...) numa nova aba. */
  public boolean openFile( String path ) { return mySeaView.openFile( path ); }

  public void setStatus( String text ) { mySeaView.setStatusText( text ); }

  /**
   * Adiciona um item ao menu "Plugins".
   * Exceções lançadas pela ação são mostradas num diálogo em vez de travar o programa.
   * @param accelerator atalho de teclado, ex. "ctrl shift I" (ou null)
   */
  public JMenuItem addMenuItem( String label, String accelerator, Runnable action ) {
    JMenuItem item = new JMenuItem( label );
    if( accelerator != null ) item.setAccelerator( KeyStroke.getKeyStroke( accelerator ) );
    item.addActionListener( e -> runSafely( label, action ) );
    myMenu.add( item );
    return item;
  }

  public JMenuItem addMenuItem( String label, Runnable action ) {
    return addMenuItem( label, null, action );
  }

  /**
   * Como {@link #addMenuItem}, mas só executa se houver painel ativo; senão avisa o usuário.
   */
  public JMenuItem addBundleMenuItem( String label, String accelerator, Consumer<csSeisPaneBundle> action ) {
    return addMenuItem( label, accelerator, () -> {
      csSeisPaneBundle b = getActiveBundle();
      if( b == null || b.getTraceBuffer() == null || b.getTraceBuffer().numTraces() == 0 ) {
        info( "Abra um arquivo sísmico primeiro." );
        return;
      }
      action.accept( b );
    });
  }

  public void info( String message ) {
    JOptionPane.showMessageDialog( mySeaView, message, myPluginName, JOptionPane.INFORMATION_MESSAGE );
  }

  public void error( String message ) {
    JOptionPane.showMessageDialog( mySeaView, message, myPluginName, JOptionPane.ERROR_MESSAGE );
  }

  private void runSafely( String label, Runnable action ) {
    try {
      action.run();
    }
    catch( Throwable t ) {
      t.printStackTrace();
      error( label + " falhou:\n" + t );
    }
  }
}
