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
 * Non-modal window showing cseis/resources/help_clock.html (Help menu and the 'Ajuda' button of the
 * clock plugin). Only one window is kept open.
 */
public class csClockHelp {
  public static final String TITLE = "Ajuda - Clock dos nodes (OBN)";
  private static JDialog ourDialog = null;

  public static void show( Component parent ) {
    Window owner = ( parent instanceof Window ) ? (Window)parent
                 : ( parent == null ) ? null : SwingUtilities.getWindowAncestor( parent );
    if( ourDialog != null && ourDialog.isDisplayable() ) {
      if( ourDialog.getOwner() == owner ) {
        ourDialog.toFront();
        return;
      }
      ourDialog.dispose();   // reopen with the new owner (a window owned by a hidden modal dialog would be blocked)
    }
    JDialog dlg = new JDialog( owner, TITLE );
    dlg.setModalExclusionType( java.awt.Dialog.ModalExclusionType.APPLICATION_EXCLUDE );
    JEditorPane pane = new JEditorPane();
    pane.setContentType( "text/html" );
    pane.setEditable( false );
    pane.setText( loadHtml() );
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
    ourDialog = dlg;
  }

  static String loadHtml() {
    try( InputStream in = csClockHelp.class.getResourceAsStream( "/cseis/resources/help_clock.html" ) ) {
      if( in == null ) return "<html><body>help_clock.html não encontrado no SeaView.jar</body></html>";
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
