/* SeaView plugin: f-c dispersion image of the active pane (e.g. a VSG from Interferometria)
   M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicTraceBuffer;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.util.Locale;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 * Phase velocity vs frequency (phase-shift method) of the traces of the active pane, using a header
 * with the source-receiver distance (offset). Intended for virtual shot gathers (correlations), but
 * works for any gather whose time axis starts at the source time ("lag 0").
 */
public class csPluginDispersao implements csISeaViewPlugin {
  private final csDispersionFC.Params myParams = new csDispersionFC.Params();
  private int myWinChoice = 0;
  private String myOffHdr = "offset";

  @Override
  public String getName() {
    return "Dispersão f-c (painel ativo)...";
  }
  @Override
  public String getDescription() {
    return "Imagem velocidade de fase x frequência (phase-shift) do painel ativo, p.ex. um VSG da Interferometria, usando o offset";
  }

  @Override
  public void run( csIPluginContext ctx ) {
    csISeismicTraceBuffer buf = ctx.getTraceBuffer();
    csHeaderDef[] defs = ctx.getHeaderDef();
    if( buf == null || buf.numTraces() < 2 || defs == null ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "São necessários pelo menos 2 traços no painel ativo", "Dispersão f-c", JOptionPane.WARNING_MESSAGE );
      return;
    }
    double dt = ctx.getSampleInt() / 1000.0;
    int ns = buf.numSamples();
    String[] names = new String[defs.length];
    for( int i = 0; i < defs.length; i++ ) names[i] = defs[i].name;
    if( ctx.getHeaderIndex( myOffHdr ) < 0 ) myOffHdr = names[0];

    // lag 0: from the VSG title ("lag 0 em ... ms"); VSG with both sides summed starts at lag 0
    double lag0ms = csDispersionFC.lag0FromTitle( ctx.getTitle() );
    boolean summedVsg = ctx.getTitle() != null && ctx.getTitle().contains( "lados somados" );
    // VSG written by the SeaSeis module INTERFEROMETRIA: headers lag0_ms and xc_sides
    int ihLag0 = ctx.getHeaderIndex( "lag0_ms" ), ihSides = ctx.getHeaderIndex( "xc_sides" );
    if( ihLag0 >= 0 ) lag0ms = buf.headerValues( 0 )[ihLag0].doubleValue();
    if( ihSides >= 0 && buf.headerValues( 0 )[ihSides].intValue() == 1 ) summedVsg = true;
    if( lag0ms < 0 ) lag0ms = 0;
    if( myParams.fmax <= 0 || myParams.fmax > 0.5 / dt ) myParams.fmax = Math.min( 2.5, 0.5 / dt );

    JComboBox<String> comboOff = new JComboBox<String>( names );
    comboOff.setSelectedItem( myOffHdr );
    comboOff.setMaximumRowCount( 20 );
    JTextField tLag0 = new JTextField( csDispersionFC.fmt( lag0ms / 1000.0 ), 7 );
    JComboBox<String> comboSide = new JComboBox<String>( csDispersionFC.SIDES );
    comboSide.setSelectedIndex( myParams.side );
    JTextField tTmax = new JTextField( csDispersionFC.fmt( myParams.tmax ), 6 );
    JComboBox<String> comboWin = new JComboBox<String>( csDispersionFC.WINDOWS );
    comboWin.setSelectedIndex( myWinChoice );
    JTextField tVcut = new JTextField( csDispersionFC.fmt( myParams.vcut ), 6 ), tPad = new JTextField( csDispersionFC.fmt( myParams.pad ), 4 );
    comboWin.setToolTipText( "Teste de alias: um ramo real rápido aparece na janela rápida; um alias de onda lenta aparece com velocidade alta mas a energia está na janela lenta" );
    tVcut.setToolTipText( "Velocidade de corte [m/s]: entre a velocidade de grupo das ondas rápidas e a das lentas (linha 1281: modos da água ~1270, Scholte ~470 -> 1200)" );
    tPad.setToolTipText( "Folga [s] para a duração do pulso; a fronteira é suavizada em +-0,5 s" );
    JTextField tFmin = new JTextField( csDispersionFC.fmt( myParams.fmin ), 6 ), tFmax = new JTextField( csDispersionFC.fmt( myParams.fmax ), 6 );
    JTextField tCmin = new JTextField( csDispersionFC.fmt( myParams.cmin ), 6 ), tCmax = new JTextField( csDispersionFC.fmt( myParams.cmax ), 6 );
    JTextField tDmin = new JTextField( csDispersionFC.fmt( myParams.dmin ), 6 );
    JTextField tDmax = new JTextField( myParams.dmax >= 1e29 ? "" : csDispersionFC.fmt( myParams.dmax ), 6 );
    tLag0.setToolTipText( "Tempo (eixo vertical) do lag 0. VSG com os dois lados: o título do painel informa (lag 0 em ... ms). Lados somados: 0" );
    tTmax.setToolTipText( "Maior lag usado [s]. 0 = até o fim do traço. O offset máximo útil é ~ c_min x lag máximo" );
    tDmax.setToolTipText( "Vazio = sem limite" );
    comboSide.setEnabled( !summedVsg );

    JPanel p = new JPanel( new GridBagLayout() );
    int row = 0;
    row = addRow( p, row, "Cabeçalho de distância", comboOff, new JLabel( " (valor absoluto)" ), null );
    row = addRow( p, row, "Lag 0 no tempo [s]", tLag0, new JLabel( " lado" ), comboSide );
    row = addRow( p, row, "Lag máximo usado [s]", tTmax, new JLabel( " (0 = traço todo)" ), null );
    row = addRow( p, row, "Frequência [Hz]  mín", tFmin, new JLabel( " máx" ), tFmax );
    row = addRow( p, row, "Velocidade [m/s]  mín", tCmin, new JLabel( " máx" ), tCmax );
    row = addRow( p, row, "Offset [m]  mín", tDmin, new JLabel( " máx" ), tDmax );
    row = addRow( p, row, "Janela no tempo", comboWin, null, null );
    row = addRow( p, row, "   corte: v [m/s]", tVcut, new JLabel( " folga [s]" ), tPad );
    if( summedVsg ) row = addRow( p, row, "", new JLabel( "VSG com os lados já somados: usa o traço a partir do lag 0" ), null, null );
    if( JOptionPane.showConfirmDialog( ctx.getParentFrame(), p, "Dispersão f-c - " + ctx.getTitle(), JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE ) != JOptionPane.OK_OPTION ) return;

    double lag0;
    try {
      lag0 = num( tLag0 );
      myParams.tmax = num( tTmax );
      myParams.vcut = num( tVcut );
      myParams.pad  = num( tPad );
      myParams.fmin = num( tFmin ); myParams.fmax = num( tFmax );
      myParams.cmin = num( tCmin ); myParams.cmax = num( tCmax );
      myParams.dmin = num( tDmin );
      myParams.dmax = tDmax.getText().trim().isEmpty() ? 1.0e30 : num( tDmax );
    }
    catch( NumberFormatException e ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "Número inválido", "Dispersão f-c", JOptionPane.ERROR_MESSAGE );
      return;
    }
    myOffHdr = (String)comboOff.getSelectedItem();
    myParams.side = comboSide.getSelectedIndex();
    myParams.oneSided = summedVsg;
    myParams.lag0 = (int)Math.round( lag0 / dt );
    // lag 0 on the first sample: there are no negative lags, so the trace is one-sided (e.g. a VSG with
    // both sides already summed that has no lag0_ms/xc_sides headers and whose title is not the original one)
    if( !myParams.oneSided && myParams.lag0 == 0 && myParams.side != csDispersionFC.SIDE_CAUSAL ) myParams.oneSided = true;
    if( myParams.lag0 < 0 || myParams.lag0 >= ns ) { error( ctx, "Lag 0 fora do traço" ); return; }
    if( myParams.cmax <= myParams.cmin || myParams.cmin <= 0 ) { error( ctx, "Velocidades: 0 < mín < máx" ); return; }
    if( myParams.fmax <= myParams.fmin ) { error( ctx, "fmax deve ser > fmin" ); return; }

    int ih = ctx.getHeaderIndex( myOffHdr );
    int ntr = buf.numTraces();
    float[][] tr = new float[ntr][];
    double[] offs = new double[ntr];
    for( int i = 0; i < ntr; i++ ) {
      tr[i] = buf.samples( i );
      offs[i] = buf.headerValues( i )[ih].doubleValue();
    }
    myWinChoice = comboWin.getSelectedIndex();
    if( myWinChoice != 0 && myParams.vcut <= 0 ) { error( ctx, "A velocidade de corte deve ser > 0" ); return; }
    if( myWinChoice == 3 ) {
      for( int w : new int[]{ csDispersionFC.WIN_FAST, csDispersionFC.WIN_SLOW } ) {
        csDispersionFC.Params q = copyParams( myParams );
        q.timeWindow = w;
        show( ctx, tr, offs, dt, q, ctx.getTitle(), myOffHdr );
      }
    }
    else {
      myParams.timeWindow = myWinChoice;
      show( ctx, tr, offs, dt, myParams, ctx.getTitle(), myOffHdr );
    }
  }

  static csDispersionFC.Params copyParams( csDispersionFC.Params a ) {
    csDispersionFC.Params b = new csDispersionFC.Params();
    b.lag0 = a.lag0; b.side = a.side; b.oneSided = a.oneSided; b.tmax = a.tmax; b.fmin = a.fmin; b.fmax = a.fmax;
    b.cmin = a.cmin; b.cmax = a.cmax; b.nc = a.nc; b.dmin = a.dmin; b.dmax = a.dmax;
    b.timeWindow = a.timeWindow; b.vcut = a.vcut; b.pad = a.pad; b.edge = a.edge;
    return b;
  }

  /** Computes and opens the f-c window (also used by Interferometria) */
  static void show( csIPluginContext ctx, float[][] traces, double[] offs, double dt, csDispersionFC.Params prm, String title, String offName ) {
    if( prm.timeWindow == csDispersionFC.WIN_FAST ) title = title + String.format( Locale.US, " [janela rápida: t < x/%.0f + %.3g s]", prm.vcut, prm.pad );
    else if( prm.timeWindow == csDispersionFC.WIN_SLOW ) title = title + String.format( Locale.US, " [janela lenta: t > x/%.0f + %.3g s]", prm.vcut, prm.pad );
    csDispersionFC.Image img;
    try {
      img = csDispersionFC.compute( traces, offs, dt, prm );
    }
    catch( IllegalArgumentException e ) {
      error( ctx, e.getMessage() );
      return;
    }
    double tUsed = ( img.nLagSamples - 1 ) * dt;
    String side = prm.oneSided ? "traço a partir do lag 0" : csDispersionFC.SIDES[prm.side];
    StringBuilder info = new StringBuilder();
    info.append( String.format( Locale.US, "<b>%s</b><br>%d traços, %s de %.0f a %.0f m (Δx típico %.0f m), %s, lags 0-%.4g s, Δf = %.4g Hz",
        title, img.nUsed, offName, img.dmin, img.dmax, img.dx, side, tUsed, img.freqs.length > 1 ? img.freqs[1] - img.freqs[0] : 0.0 ) );
    // a wave slower than dmax / tmax does not reach the farthest traces inside the lag window
    double cTrunc = img.dmax / Math.max( tUsed, dt );
    if( cTrunc > prm.cmin ) {
      info.append( String.format( Locale.US, "<br><font color='#b00000'>Atenção: com lags até %.4g s, ondas mais lentas que %.0f m/s não chegam ao offset %.0f m dentro da janela; reduza o offset máximo ou aumente o lag</font>",
          tUsed, cTrunc, img.dmax ) );
    }
    csDispersionFC.Frame f = new csDispersionFC.Frame( img, title, info.toString() );
    f.setLocationRelativeTo( ctx.getParentFrame() );
    if( prm.timeWindow == csDispersionFC.WIN_SLOW ) f.setLocation( f.getX() + 60, f.getY() + 60 );
    f.setVisible( true );
  }

  private static double num( JTextField t ) { return Double.parseDouble( t.getText().trim().replace( ',', '.' ) ); }
  private static void error( csIPluginContext ctx, String text ) {
    JOptionPane.showMessageDialog( ctx.getParentFrame(), text, "Dispersão f-c", JOptionPane.ERROR_MESSAGE );
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
