/* Help window for the OBN clock plugins. M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Window;
import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JEditorPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;

/**
 * Non-modal help windows with pages of cseis/resources: help_clock.html (Help menu and the 'Ajuda' button of
 * the clock plugin) and help_difusividade.html (plugin Difusividade do campo). One window per page.
 */
public class csClockHelp {
  public static final String TITLE = "Ajuda - Clock dos nodes (OBN)";
  public static final String TITLE_DIFUSIVIDADE = "Ajuda - Difusividade do campo (SPAC, f-k, simetria, atrasos)";
  /** one open window per help page */
  private static final java.util.Map<String,JDialog> ourDialogs = new java.util.HashMap<String,JDialog>();

  public static void show( Component parent ) {
    show( parent, TITLE, "help_clock.html" );
  }
  /** Help of the plugin Difusividade do campo */
  public static void showDifusividade( Component parent ) {
    show( parent, TITLE_DIFUSIVIDADE, "help_difusividade.html" );
  }

  /** Non-modal window with cseis/resources/'resource' */
  public static void show( Component parent, String title, String resource ) {
    Window owner = ( parent instanceof Window ) ? (Window)parent
                 : ( parent == null ) ? null : SwingUtilities.getWindowAncestor( parent );
    JDialog old = ourDialogs.get( resource );
    if( old != null && old.isDisplayable() ) {
      if( old.getOwner() == owner ) {
        old.toFront();
        return;
      }
      old.dispose();   // reopen with the new owner (a window owned by a hidden modal dialog would be blocked)
    }
    JDialog dlg = new JDialog( owner, title );
    dlg.setModalExclusionType( java.awt.Dialog.ModalExclusionType.APPLICATION_EXCLUDE );
    JEditorPane pane = new JEditorPane();
    pane.setContentType( "text/html" );
    pane.setEditable( false );
    pane.setText( loadHtml( resource ) );
    pane.setCaretPosition( 0 );
    JScrollPane scroll = new JScrollPane( pane );
    scroll.setPreferredSize( new Dimension( 820, 700 ) );
    JButton close = new JButton( "Fechar" );
    close.addActionListener( e -> dlg.dispose() );
    JPanel south = new JPanel();
    south.add( close );
    dlg.getContentPane().add( scroll, BorderLayout.CENTER );
    dlg.getContentPane().add( south, BorderLayout.SOUTH );
    dlg.setDefaultCloseOperation( JDialog.DISPOSE_ON_CLOSE );
    dlg.pack();
    dlg.setLocationRelativeTo( parent );
    dlg.setVisible( true );
    ourDialogs.put( resource, dlg );
  }

  static String loadHtml( String resource ) {
    try( InputStream in = csClockHelp.class.getResourceAsStream( "/cseis/resources/" + resource ) ) {
      if( in == null ) return "<html><body>" + resource + " não encontrado no SeaView.jar</body></html>";
      ByteArrayOutputStream out = new ByteArrayOutputStream();
      byte[] buf = new byte[8192];
      int n;
      while( ( n = in.read( buf ) ) > 0 ) out.write( buf, 0, n );
      return new String( out.toByteArray(), StandardCharsets.UTF_8 );
    }
    catch( Exception e ) {
      return "<html><body>Erro ao ler a ajuda: " + e + "</body></html>";
    }
  }
}
