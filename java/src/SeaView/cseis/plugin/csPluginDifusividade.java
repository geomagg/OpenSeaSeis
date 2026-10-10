/* SeaView plugin: diffuse-field tests of the noise in the active pane (SPAC / Bessel J0, f-k, hourly symmetry, neighbour delays)
   M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicTraceBuffer;
import java.awt.Cursor;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 * Raw noise of a line of receivers (one trace per receiver, long records, as used by Interferometria):
 * coherency between every pair of receivers, accumulated per time block, and four tests of the diffuse field.
 * Needs rec_x / rec_y. All traces must start at the same time.
 */
public class csPluginDifusividade implements csISeaViewPlugin {
  static final String NONE = "(todos)";
  private final csDifusividade.Params myParams = new csDifusividade.Params();
  private String mySelHdr = null, mySelVal = null, myLabelHdr = null;

  @Override
  public String getName() {
    return "Difusividade do campo (Bessel, f-k, simetria por hora)...";
  }
  @Override
  public String getDescription() {
    return "Testes de campo difuso no ruído bruto do painel ativo: coerência x distância (SPAC, Bessel J0), f-k, simetria hora a hora e atraso entre receptores vizinhos (relógio)";
  }

  @Override
  public void run( csIPluginContext ctx ) {
    csISeismicTraceBuffer buf = ctx.getTraceBuffer();
    csHeaderDef[] defs = ctx.getHeaderDef();
    if( buf == null || buf.numTraces() < 3 || defs == null ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "São necessários pelo menos 3 traços (receptores) no painel ativo", "Difusividade", JOptionPane.WARNING_MESSAGE );
      return;
    }
    int ihX = ctx.getHeaderIndex( "rec_x" ), ihY = ctx.getHeaderIndex( "rec_y" );
    if( ihX < 0 || ihY < 0 ) {
      error( ctx, "O painel precisa dos cabeçalhos rec_x e rec_y (coordenadas dos receptores).\nNo npy2cseis.py use --xy header_offset.txt." );
      return;
    }
    double dt = ctx.getSampleInt() / 1000.0;
    int ns = buf.numSamples();
    double fnyq = 0.5 / dt;
    String[] names = new String[defs.length];
    for( int i = 0; i < defs.length; i++ ) names[i] = defs[i].name;
    String[] names2 = new String[defs.length + 1];
    names2[0] = NONE;
    System.arraycopy( names, 0, names2, 1, names.length );

    if( mySelHdr == null || ( !NONE.equals( mySelHdr ) && ctx.getHeaderIndex( mySelHdr ) < 0 ) ) {
      mySelHdr = ctx.getHeaderIndex( "chan" ) >= 0 ? "chan" : NONE;
      mySelVal = null;
    }
    if( mySelVal == null && !NONE.equals( mySelHdr ) ) mySelVal = buf.headerValues( 0 )[ctx.getHeaderIndex( mySelHdr )].toString();
    if( myLabelHdr == null || ctx.getHeaderIndex( myLabelHdr ) < 0 ) {
      myLabelHdr = ctx.getHeaderIndex( "rcv" ) >= 0 ? "rcv" : ctx.getHeaderIndex( "trcno" ) >= 0 ? "trcno" : names[0];
    }
    csDifusividade.Params p = myParams;
    if( p.fmax > fnyq ) p.fmax = fnyq;
    if( p.winLen > ns * dt ) p.winLen = Math.floor( ns * dt );

    JComboBox<String> cSel = new JComboBox<String>( names2 );
    cSel.setSelectedItem( mySelHdr );
    cSel.setMaximumRowCount( 20 );
    JTextField tSelVal = new JTextField( mySelVal == null ? "" : mySelVal, 6 );
    cSel.addActionListener( e -> {
      int ih = ctx.getHeaderIndex( (String)cSel.getSelectedItem() );
      tSelVal.setText( ih >= 0 ? buf.headerValues( 0 )[ih].toString() : "" );
    } );
    JComboBox<String> cLab = new JComboBox<String>( names );
    cLab.setSelectedItem( myLabelHdr );
    cLab.setMaximumRowCount( 20 );
    JComboBox<String> cNorm = new JComboBox<String>( csPluginInterferometria.NORMS );
    cNorm.setSelectedIndex( p.tempNorm );
    JTextField tRam = tf( p.ramWin ), tWin = tf( p.winLen ), tBlk = tf( p.blockLen );
    JTextField tFmin = tf( p.fmin ), tFmax = tf( p.fmax ), tF1 = tf( p.f1 ), tF2 = tf( p.f2 ), tCmin = tf( p.cmin ), tCmax = tf( p.cmax );
    JTextField tFd1 = tf( p.fd1 ), tFd2 = tf( p.fd2 ), tCf1 = tf( p.cFitMin ), tCf2 = tf( p.cFitMax );
    JCheckBox bWhite = new JCheckBox( "Branqueamento espectral", p.whiten );
    JCheckBox bFree = new JCheckBox( "Ajustar a·J₀ + b (absorve energia vertical e ruído incoerente)", p.freeAmp );
    tWin.setToolTipText( "Janela de cálculo do espectro cruzado [s] (como as janelas de empilhamento da Interferometria)" );
    tBlk.setToolTipText( "Duração de cada bloco para a simetria e os atrasos [s] (3600 = hora a hora; 0 = registro todo)" );
    tF1.setToolTipText( "Banda da simetria f-k (linha 1281, Scholte: 0,15-0,45 Hz)" );
    tCmin.setToolTipText( "Velocidades da janela da simetria no f-k (Scholte: 450-1600 m/s; modos da água: 1500-3000 m/s)" );
    tFd1.setToolTipText( "Banda dos atrasos entre vizinhos: onde a coerência entre receptores próximos é alta (linha 1281: 0,10-0,35 Hz)" );
    bFree.setToolTipText( "Desligado: compara J0 e cos puros (teste de difusividade). Ligado: melhor para estimar c quando a coerência fica positiva longe (energia chegando na vertical)" );

    JPanel pan = new JPanel( new GridBagLayout() );
    int row = 0;
    row = addRow( pan, row, "Usar só traços com", cSel, new JLabel( " = " ), tSelVal );
    row = addRow( pan, row, "Rótulo dos receptores", cLab, null, null );
    row = addRow( pan, row, "Normalização temporal", cNorm, new JLabel( " janela [s]" ), tRam );
    row = addRow( pan, row, "Banda de pré-processamento [Hz]", tFmin, new JLabel( " até" ), tFmax );
    row = addRow( pan, row, "", bWhite, null, null );
    row = addRow( pan, row, "Janela [s]", tWin, new JLabel( " bloco [s]" ), tBlk );
    row = addRow( pan, row, "Simetria: banda [Hz]", tF1, new JLabel( " até" ), tF2 );
    row = addRow( pan, row, "Simetria: velocidade [m/s]", tCmin, new JLabel( " até" ), tCmax );
    row = addRow( pan, row, "Atrasos entre vizinhos: banda [Hz]", tFd1, new JLabel( " até" ), tFd2 );
    row = addRow( pan, row, "SPAC: busca de c [m/s]", tCf1, new JLabel( " até" ), tCf2 );
    row = addRow( pan, row, "", bFree, null, null );
    JPanel note = new JPanel( new FlowLayout( FlowLayout.LEFT, 0, 0 ) );
    note.add( new JLabel( "<html><i>Painel com o ruído bruto: um traço por receptor, todos começando no mesmo instante.</i></html>" ) );
    row = addRow( pan, row, "", note, null, null );
    if( JOptionPane.showConfirmDialog( ctx.getParentFrame(), pan, "Difusividade do campo - " + ctx.getTitle(), JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE ) != JOptionPane.OK_OPTION ) return;

    try {
      p.ramWin = num( tRam ); p.winLen = num( tWin ); p.blockLen = num( tBlk );
      p.fmin = num( tFmin ); p.fmax = num( tFmax ); p.f1 = num( tF1 ); p.f2 = num( tF2 );
      p.cmin = num( tCmin ); p.cmax = num( tCmax ); p.fd1 = num( tFd1 ); p.fd2 = num( tFd2 );
      p.cFitMin = num( tCf1 ); p.cFitMax = num( tCf2 );
    }
    catch( NumberFormatException e ) {
      error( ctx, "Número inválido" );
      return;
    }
    p.tempNorm = cNorm.getSelectedIndex();
    p.whiten = bWhite.isSelected();
    p.freeAmp = bFree.isSelected();
    mySelHdr = (String)cSel.getSelectedItem();
    mySelVal = tSelVal.getText().trim();
    myLabelHdr = (String)cLab.getSelectedItem();
    if( p.fmax <= p.fmin ) { error( ctx, "Banda de pré-processamento: mín < máx" ); return; }
    if( p.cmin <= 0 || p.cmax <= p.cmin || p.cFitMin <= 0 || p.cFitMax <= p.cFitMin ) { error( ctx, "Velocidades: 0 < mín < máx" ); return; }
    if( p.winLen * 4 > ns * dt ) {
      int ok = JOptionPane.showConfirmDialog( ctx.getParentFrame(), String.format( Locale.US, "O traço tem só %.4g janelas de %.4g s. A coerência fica pouco estável. Continuar?", ns * dt / p.winLen, p.winLen ),
          "Difusividade", JOptionPane.OK_CANCEL_OPTION, JOptionPane.WARNING_MESSAGE );
      if( ok != JOptionPane.OK_OPTION ) return;
    }

    //--- traces
    int ihSel = NONE.equals( mySelHdr ) ? -1 : ctx.getHeaderIndex( mySelHdr );
    int ihLab = ctx.getHeaderIndex( myLabelHdr );
    int ihT1 = ctx.getHeaderIndex( "time_samp1" ), ihT1us = ctx.getHeaderIndex( "time_samp1_us" );
    List<Integer> use = new ArrayList<Integer>();
    for( int i = 0; i < buf.numTraces(); i++ ) {
      if( ihSel >= 0 && !matches( buf.headerValues( i )[ihSel].toString(), mySelVal ) ) continue;
      use.add( i );
    }
    if( use.size() < 3 ) { error( ctx, "Menos de 3 traços com " + mySelHdr + " = " + mySelVal ); return; }
    int n = use.size();
    float[][] tr = new float[n][];
    double[] x = new double[n], y = new double[n];
    String[] lab = new String[n];
    double tmin = Double.POSITIVE_INFINITY, tmax = Double.NEGATIVE_INFINITY;
    for( int k = 0; k < n; k++ ) {
      int i = use.get( k );
      tr[k] = buf.samples( i );
      x[k] = buf.headerValues( i )[ihX].doubleValue();
      y[k] = buf.headerValues( i )[ihY].doubleValue();
      lab[k] = buf.headerValues( i )[ihLab].toString().trim();
      if( ihT1 >= 0 ) {
        double t = buf.headerValues( i )[ihT1].doubleValue() + ( ihT1us >= 0 ? buf.headerValues( i )[ihT1us].doubleValue() * 1e-6 : 0 );
        tmin = Math.min( tmin, t ); tmax = Math.max( tmax, t );
      }
    }

    csDifusividade.Result res;
    ctx.getParentFrame().setCursor( Cursor.getPredefinedCursor( Cursor.WAIT_CURSOR ) );
    try {
      res = csDifusividade.compute( tr, x, y, lab, dt, p );
    }
    catch( IllegalArgumentException e ) {
      error( ctx, e.getMessage() );
      return;
    }
    finally {
      ctx.getParentFrame().setCursor( Cursor.getDefaultCursor() );
    }
    if( ihT1 >= 0 && tmax - tmin > 0.5 * dt ) {
      res.warnings += String.format( Locale.US, "Os traços não começam no mesmo instante (time_samp1 varia %.4g s): alinhe antes (ALIGN_TIME), senão os atrasos e a coerência ficam errados. ", tmax - tmin );
    }
    String title = ctx.getTitle() + ( ihSel >= 0 ? " (" + mySelHdr + " " + mySelVal + ")" : "" );
    csDifusividadeFrame f = new csDifusividadeFrame( res, title );
    f.setLocationRelativeTo( ctx.getParentFrame() );
    f.setVisible( true );
  }

  private static boolean matches( String h, String v ) {
    try { return Math.abs( Double.parseDouble( h.trim() ) - Double.parseDouble( v.trim() ) ) < 1e-6; }
    catch( NumberFormatException e ) { return h.trim().equals( v.trim() ); }
  }
  private static JTextField tf( double v ) { return new JTextField( csDispersionFC.fmt( v ), 6 ); }
  private static double num( JTextField t ) { return Double.parseDouble( t.getText().trim().replace( ',', '.' ) ); }
  private static void error( csIPluginContext ctx, String text ) {
    JOptionPane.showMessageDialog( ctx.getParentFrame(), text, "Difusividade", JOptionPane.ERROR_MESSAGE );
  }
  private static int addRow( JPanel p, int row, String label, JComponent c1, JComponent c2, JComponent c3 ) {
    GridBagConstraints g = new GridBagConstraints();
    g.insets = new Insets( 3, 4, 3, 4 );
    g.anchor = GridBagConstraints.WEST;
    g.gridy = row;
    g.gridx = 0; p.add( new JLabel( label ), g );
    g.gridx = 1; if( c1 != null ) p.add( c1, g );
    g.gridx = 2; if( c2 != null ) p.add( c2, g );
    g.gridx = 3; if( c3 != null ) p.add( c3, g );
    return row + 1;
  }
}
