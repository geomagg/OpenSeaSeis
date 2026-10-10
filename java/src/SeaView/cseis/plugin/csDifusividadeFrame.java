/* Window of the diffuse-field tests (csDifusividade): SPAC / Bessel J0, f-k, hourly symmetry, neighbour delays
   M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;

public final class csDifusividadeFrame extends JFrame {
  static final Color GREEN = new Color( 47, 125, 60 ), RED = new Color( 190, 50, 40 ), ORANGE = new Color( 205, 120, 20 ),
                     BLUE = new Color( 15, 95, 122 ), GREY = new Color( 165, 165, 165 ), MAGENTA = new Color( 220, 60, 200 );

  private final csDifusividade.Result myR;
  private final JTabbedPane myTabs = new JTabbedPane();
  private final JLabel myStatus = new JLabel( " " );
  // SPAC tab state
  private int myG;
  private double myC;
  private final JComboBox<String> myFreqBox;
  private final JTextField myCText = new JTextField( 6 );
  private final JLabel mySpacInfo = new JLabel( " " );
  private Chart myCohChart;

  public csDifusividadeFrame( csDifusividade.Result r, String title ) {
    super( "Difusividade do campo — " + title );
    myR = r;
    setDefaultCloseOperation( DISPOSE_ON_CLOSE );
    String[] fs = new String[r.freqs.length];
    for( int g = 0; g < fs.length; g++ ) fs[g] = String.format( Locale.US, "%.3f Hz", r.freqs[g] );
    myFreqBox = new JComboBox<String>( fs );
    myFreqBox.setMaximumRowCount( 25 );
    myG = nearest( r.freqs, 0.5 * ( r.p.f1 + r.p.f2 ) );
    myC = r.cFit[myG];

    JLabel head = new JLabel( header() );
    head.setBorder( BorderFactory.createEmptyBorder( 6, 10, 4, 10 ) );
    JButton csv = new JButton( "Exportar CSVs..." ), png = new JButton( "Salvar PNG da aba..." ), help = new JButton( "Ajuda" );
    help.setToolTipText( "Como ler cada aba (também no menu Help)" );
    help.addActionListener( e -> csClockHelp.showDifusividade( this ) );
    csv.addActionListener( e -> exportCsv() );
    png.addActionListener( e -> savePng() );
    JPanel bar = new JPanel( new FlowLayout( FlowLayout.LEFT, 8, 2 ) );
    bar.add( help ); bar.add( csv ); bar.add( png ); bar.add( myStatus );
    JPanel top = new JPanel( new BorderLayout() );
    top.add( head, BorderLayout.CENTER ); top.add( bar, BorderLayout.SOUTH );

    myTabs.addTab( "Bessel J₀ (SPAC)", spacTab() );
    myTabs.addTab( "f–k", fkTab() );
    myTabs.addTab( "Simetria por hora", blockTab() );
    myTabs.addTab( "Atraso entre vizinhos (relógio)", delayTab() );
    getContentPane().setLayout( new BorderLayout() );
    getContentPane().add( top, BorderLayout.NORTH );
    getContentPane().add( myTabs, BorderLayout.CENTER );
    setPreferredSize( new Dimension( 1060, 820 ) );
    pack();
    selectFreq( myG, Double.NaN );
  }

  private String header() {
    csDifusividade.Result r = myR;
    StringBuilder sb = new StringBuilder( "<html>" );
    sb.append( String.format( Locale.US, "<b>%d receptores</b>, Δx típico %.0f m, linha no azimute %.0f° (+s do 1º ao último traço) &nbsp;·&nbsp; %d janelas de %s s em %d bloco(s) de %s s",
        r.n, r.dx, r.azimuth, r.nwinTotal, csDispersionFC.fmt( r.p.winLen ), r.nblk, csDispersionFC.fmt( r.p.blockLen > 0 ? r.p.blockLen : r.nwinTotal * r.p.winLen ) ) );
    sb.append( String.format( Locale.US, "<br>Pré-processamento %s–%s Hz%s%s &nbsp;·&nbsp; banda de análise %s–%s Hz, %.0f–%.0f m/s &nbsp;·&nbsp; atrasos %s–%s Hz",
        csDispersionFC.fmt( r.p.fmin ), csDispersionFC.fmt( r.p.fmax ), r.p.tempNorm == 2 ? String.format( Locale.US, ", média absoluta móvel %.3g s", r.p.ramWin ) : r.p.tempNorm == 1 ? ", one-bit" : "",
        r.p.whiten ? ", branqueado" : "", csDispersionFC.fmt( r.p.f1 ), csDispersionFC.fmt( r.p.f2 ), r.p.cmin, r.p.cmax,
        csDispersionFC.fmt( r.p.fd1 ), csDispersionFC.fmt( r.p.fd2 ) ) );
    if( !r.warnings.isEmpty() ) sb.append( "<br><font color='#b00000'>" ).append( r.warnings ).append( "</font>" );
    return sb.append( "</html>" ).toString();
  }

  static int nearest( double[] v, double x ) {
    int b = 0;
    for( int i = 1; i < v.length; i++ ) if( Math.abs( v[i] - x ) < Math.abs( v[b] - x ) ) b = i;
    return b;
  }
  static String rumo( double az ) { return csAssimetria.rumo( ( az % 360 + 360 ) % 360 ); }

  //====================================================================
  // 1. SPAC
  //====================================================================
  private JComponent spacTab() {
    csDifusividade.Result r = myR;
    myCohChart = new Chart( "distância entre receptores r (km)", "coerência", myStatus );
    myCohChart.setRange( 0, maxDist() / 1000.0 * 1.02, -1, 1 );
    myCohChart.hLines.add( 0.0 );

    // misfit image, normalized per frequency (1 = best c)
    int nf = r.freqs.length, nc = r.cGrid.length;
    double[][] v = new double[nf][nc];
    for( int g = 0; g < nf; g++ ) {
      double lo = Double.POSITIVE_INFINITY, hi = Double.NEGATIVE_INFINITY;
      for( double e : r.spacMisfit[g] ) if( !Double.isNaN( e ) ) { lo = Math.min( lo, e ); hi = Math.max( hi, e ); }
      for( int k = 0; k < nc; k++ ) {
        double e = r.spacMisfit[g][k];
        v[g][k] = Double.isNaN( e ) || !( hi > lo ) ? Double.NaN : ( hi - e ) / ( hi - lo );
      }
    }
    Heat img = new Heat( r.freqs, r.cGrid, v, 0, 1, "frequência (Hz)", "velocidade de fase c (m/s)", myStatus );
    img.reader = ( f, c ) -> {
      int g = nearest( r.freqs, f ), k = nearest( r.cGrid, c );
      return String.format( Locale.US, "f = %.3f Hz, c = %.0f m/s: rms para J₀ = %.3f  (clique: escolhe f e c)", r.freqs[g], r.cGrid[k], r.spacMisfit[g][k] );
    };
    img.click = ( f, c ) -> selectFreq( nearest( r.freqs, f ), c );
    // overlays: best c (filled = J0 explains clearly better than cos), first zero, spacing limit
    List<double[]> good = new ArrayList<double[]>(), weak = new ArrayList<double[]>(), zero = new ArrayList<double[]>();
    double[] lim = new double[nf];
    for( int g = 0; g < nf; g++ ) {
      lim[g] = csDifusividade.cLow( r.freqs[g], r.dx );
      if( !Double.isNaN( r.cFit[g] ) ) ( clear( g ) ? good : weak ).add( new double[]{ r.freqs[g], r.cFit[g] } );
      if( !Double.isNaN( r.cZero[g] ) ) zero.add( new double[]{ r.freqs[g], r.cZero[g] } );
    }
    img.series.add( new Ser( r.freqs, lim, Color.white, Ser.DASH, "limite do espaçamento (1º zero de J₀ = Δx)" ) );
    img.series.add( Ser.of( good, Color.white, Ser.DOT, "melhor c (J₀ ≫ cos)" ) );
    img.series.add( Ser.of( weak, Color.white, Ser.OPEN, "melhor c (fraco)" ) );
    img.series.add( Ser.of( zero, MAGENTA, Ser.CROSS, "c pelo 1º zero" ) );

    myFreqBox.addActionListener( e -> { if( myFreqBox.getSelectedIndex() != myG ) selectFreq( myFreqBox.getSelectedIndex(), Double.NaN ); } );
    myCText.addActionListener( e -> {
      try { selectFreq( myG, Double.parseDouble( myCText.getText().trim().replace( ',', '.' ) ) ); }
      catch( NumberFormatException ex ) { myStatus.setText( "Velocidade inválida" ); }
    } );
    JButton auto = new JButton( "c ajustado" );
    auto.addActionListener( e -> selectFreq( myG, Double.NaN ) );
    JPanel ctl = new JPanel( new FlowLayout( FlowLayout.LEFT, 6, 2 ) );
    ctl.add( new JLabel( "Frequência" ) ); ctl.add( myFreqBox );
    ctl.add( new JLabel( "  c [m/s]" ) ); ctl.add( myCText ); ctl.add( auto );
    ctl.add( new JLabel( "  (Enter aplica; clique na imagem abaixo para escolher f e c)" ) );
    JPanel up = new JPanel( new BorderLayout() );
    up.add( ctl, BorderLayout.NORTH ); up.add( myCohChart, BorderLayout.CENTER ); up.add( mySpacInfo, BorderLayout.SOUTH );
    mySpacInfo.setBorder( BorderFactory.createEmptyBorder( 4, 10, 4, 10 ) );

    JLabel note = new JLabel( "<html>Imagem: ajuste da coerência a J₀(2πfr/c) para cada (f, c), normalizado por frequência (amarelo = melhor c). "
        + "Abaixo da linha tracejada o 1º zero de J₀ fica mais perto que um espaçamento e o SPAC não resolve c (alias espacial). "
        + "Pontos cheios: J₀ explica os dados muito melhor que uma onda plana (rms J₀ &lt; 0,25 e rms cos − rms J₀ &gt; 0,2).</html>" );
    note.setBorder( BorderFactory.createEmptyBorder( 2, 10, 2, 10 ) );
    JPanel down = new JPanel( new BorderLayout() );
    down.add( img, BorderLayout.CENTER ); down.add( note, BorderLayout.SOUTH );
    JSplitPane sp = new JSplitPane( JSplitPane.VERTICAL_SPLIT, up, down );
    sp.setResizeWeight( 0.5 );
    sp.setDividerLocation( 330 );
    return sp;
  }

  private boolean clear( int g ) {
    return myR.rmsJ0[g] < 0.25 && myR.rmsCos[g] - myR.rmsJ0[g] > 0.2;
  }

  private double maxDist() {
    double m = 0;
    for( double[] d : myR.dist ) for( double v : d ) m = Math.max( m, v );
    return m;
  }

  /** Selects the frequency group g and the velocity c (NaN = fitted) for the coherency chart */
  private void selectFreq( int g, double c ) {
    csDifusividade.Result r = myR;
    myG = g;
    if( myFreqBox.getSelectedIndex() != g ) myFreqBox.setSelectedIndex( g );
    boolean fitted = Double.isNaN( c );
    myC = fitted ? r.cFit[g] : c;
    myCText.setText( Double.isNaN( myC ) ? "" : String.format( Locale.US, "%.0f", myC ) );
    double f = r.freqs[g];
    int n = r.n;
    List<double[]> pre = new ArrayList<double[]>();
    double[] co = new double[2];
    for( int i = 0; i < n; i++ ) for( int j = i + 1; j < n; j++ ) {
      csDifusividade.coh( r.autoT[g], r.crossT[g], i, j, n, co );
      if( !Double.isNaN( co[0] ) ) pre.add( new double[]{ r.dist[i][j] / 1000.0, co[0] } );
    }
    double[][] bins = csDifusividade.binned( r, r.autoT[g], r.crossT[g] );
    double[] bx = new double[bins[0].length];
    for( int b = 0; b < bx.length; b++ ) bx[b] = bins[0][b] / 1000.0;
    myCohChart.series.clear();
    myCohChart.series.add( Ser.of( pre, GREY, Ser.DOT_SMALL, "Re γ, pares" ) );
    myCohChart.series.add( new Ser( bx, bins[2], ORANGE, Ser.OPEN, "Im γ · sinal(Δs), média por distância" ) );
    myCohChart.series.add( new Ser( bx, bins[1], Color.black, Ser.DOT, "Re γ, média por distância" ) );
    String txt;
    if( !Double.isNaN( myC ) ) {
      int m = 300;
      double rm = maxDist() * 1.02;
      double[] xs = new double[m], yj = new double[m], yc = new double[m];
      double[] ab = new double[]{ 1, 0 };
      if( r.p.freeAmp ) csDifusividade.misfit( bins, f, myC, 0, true, ab );
      for( int k = 0; k < m; k++ ) {
        double rr = rm * k / ( m - 1 ), kr = 2 * Math.PI * f * rr / myC;
        xs[k] = rr / 1000.0; yj[k] = ab[0] * csDifusividade.j0( kr ) + ab[1]; yc[k] = ab[0] * Math.cos( kr ) + ab[1];
      }
      myCohChart.series.add( new Ser( xs, yj, GREEN, Ser.LINE, "J₀(2πfr/c): campo difuso" ) );
      myCohChart.series.add( new Ser( xs, yc, RED, Ser.DASH, "cos(2πfr/c): onda plana ao longo da linha" ) );
      double ej = csDifusividade.misfit( bins, f, myC, 0, r.p.freeAmp, null ), ec = csDifusividade.misfit( bins, f, myC, 1, r.p.freeAmp, null );
      txt = String.format( Locale.US, "<html><b>f = %.3f Hz, c = %.0f m/s</b> (%s, λ = %.0f m) &nbsp;·&nbsp; erro rms: <b>J₀ = %.2f</b>, cos = %.2f",
          f, myC, fitted ? "ajustado" : "escolhido", myC / f, ej, ec );
      if( r.p.freeAmp ) txt += String.format( Locale.US, " &nbsp;·&nbsp; modelo a·J₀ + b: a = %.2f, b = %+.2f", ab[0], ab[1] );
      if( !Double.isNaN( r.cZero[g] ) ) txt += String.format( Locale.US, " &nbsp;·&nbsp; c pelo 1º zero = %.0f m/s", r.cZero[g] );
      txt += "<br>" + ( ec - ej > 0.2 && ej < 0.3 ? "<font color='#2f7d3c'>A coerência decai como J₀: energia chegando de muitas direções (campo difuso).</font>"
                     : ec - ej > 0.1 ? "<font color='#a8661b'>J₀ explica melhor que o cosseno, mas o ajuste é parcial.</font>"
                                     : "<font color='#9b2f2f'>J₀ não explica melhor que uma onda plana nesta frequência.</font>" );
      txt += "</html>";
    }
    else {
      txt = String.format( Locale.US, "<html><b>f = %.3f Hz</b>: sem ajuste (mínimo no limite da busca). Escolha c na imagem ou digite.</html>", f );
    }
    mySpacInfo.setText( txt );
    myCohChart.repaint();
  }

  //====================================================================
  // 2. f-k
  //====================================================================
  private JComponent fkTab() {
    csDifusividade.Result r = myR;
    int nf = r.freqs.length, nk = r.ks.length;
    double[] kkm = new double[nk];
    for( int k = 0; k < nk; k++ ) kkm[k] = r.ks[k] * 1000.0;
    double[][] v = new double[nf][nk];
    for( int g = 0; g < nf; g++ ) for( int k = 0; k < nk; k++ ) v[g][k] = r.fk[k][g];
    Heat img = new Heat( r.freqs, kkm, v, 0, 1, "frequência (Hz)", "número de onda ao longo da linha k (1/km)", myStatus );
    img.reader = ( f, k ) -> {
      int g = nearest( r.freqs, f ), ik = nearest( kkm, k );
      double va = Math.abs( kkm[ik] ) > 1e-9 ? r.freqs[g] / Math.abs( kkm[ik] ) * 1000.0 : Double.POSITIVE_INFINITY;
      return String.format( Locale.US, "f = %.3f Hz, k = %+.3f 1/km, velocidade aparente %.0f m/s, P = %.2f, energia andando para %s",
          r.freqs[g], kkm[ik], va, v[g][ik], kkm[ik] >= 0 ? String.format( Locale.US, "%.0f° (%s)", r.azimuth, rumo( r.azimuth ) )
                                                           : String.format( Locale.US, "%.0f° (%s)", ( r.azimuth + 180 ) % 360, rumo( r.azimuth + 180 ) ) );
    };
    // velocity lines of the analysis window, and the SPAC c
    for( double c : new double[]{ r.p.cmin, r.p.cmax } ) for( int sg : new int[]{ 1, -1 } ) {
      double[] y = new double[nf];
      for( int g = 0; g < nf; g++ ) { y[g] = sg * r.freqs[g] / c * 1000.0; if( Math.abs( y[g] ) > kkm[nk - 1] ) y[g] = Double.NaN; }
      img.series.add( new Ser( r.freqs, y, Color.white, Ser.DASH, sg > 0 && c == r.p.cmin ? String.format( Locale.US, "±f/%.0f e ±f/%.0f (janela da simetria)", r.p.cmin, r.p.cmax ) : null ) );
    }
    List<double[]> sp = new ArrayList<double[]>();
    for( int g = 0; g < nf; g++ ) if( !Double.isNaN( r.cFit[g] ) && clear( g ) ) {
      double k = r.freqs[g] / r.cFit[g] * 1000.0;
      if( k <= kkm[nk - 1] ) { sp.add( new double[]{ r.freqs[g], k } ); sp.add( new double[]{ r.freqs[g], -k } ); }
    }
    img.series.add( Ser.of( sp, MAGENTA, Ser.DOT_SMALL, "±f/c do SPAC" ) );
    img.vLines.add( r.p.f1 ); img.vLines.add( r.p.f2 );
    double to = r.azimuth, from = ( r.azimuth + 180 ) % 360;
    JLabel note = new JLabel( String.format( Locale.US, "<html>f–k do ruído bruto (todas as janelas), normalizado por frequência. Número de Nyquist ±%.3f 1/km (Δx = %.0f m). "
        + "<b>k &gt; 0</b>: energia andando para %.0f° (%s), isto é, chegando de %.0f° (%s); <b>k &lt; 0</b>: o contrário. "
        + "Campo difuso: energia nos dois lados, em ±f/c. Linhas verticais: banda da simetria.</html>",
        kkm[nk - 1], r.dx, to, rumo( to ), from, rumo( from ) ) );
    note.setBorder( BorderFactory.createEmptyBorder( 4, 10, 4, 10 ) );
    JPanel p = new JPanel( new BorderLayout() );
    p.add( img, BorderLayout.CENTER ); p.add( note, BorderLayout.SOUTH );
    return p;
  }

  //====================================================================
  // 3. Symmetry per block
  //====================================================================
  private JComponent blockTab() {
    csDifusividade.Result r = myR;
    int nb = r.nblk;
    double[] t = new double[nb];
    for( int b = 0; b < nb; b++ ) t[b] = ( r.blkStart[b] + 0.5 * r.blkWins[b] * r.p.winLen ) / 3600.0;
    double t1 = nb > 0 ? ( r.blkStart[nb - 1] + r.blkWins[nb - 1] * r.p.winLen ) / 3600.0 : 1;
    Chart cd = new Chart( "tempo desde o início do painel (h)", "D = (P+ − P−)/(P+ + P−)", myStatus );
    cd.setRange( 0, Math.max( t1, 1e-3 ), -1, 1 );
    cd.bandLo = -0.3; cd.bandHi = 0.3;
    cd.hLines.add( 0.0 );
    if( !Double.isNaN( r.dTotal ) ) cd.series.add( new Ser( new double[]{ 0, t1 }, new double[]{ r.dTotal, r.dTotal }, BLUE, Ser.DASH, String.format( Locale.US, "registro todo: D = %+.2f", r.dTotal ) ) );
    cd.series.add( new Ser( t, r.dBlk, BLUE, Ser.LINE_DOTS, "D por bloco" ) );
    Chart cr = new Chart( "tempo desde o início do painel (h)", "erro rms da coerência", myStatus );
    cr.setRange( 0, Math.max( t1, 1e-3 ), 0, 1 );
    cr.series.add( new Ser( t, r.rmsJ0Blk, GREEN, Ser.LINE_DOTS, "J₀ (campo difuso), c do registro todo" ) );
    cr.series.add( new Ser( t, r.rmsCosBlk, RED, Ser.LINE_DOTS, "cos (onda plana ao longo da linha), mesmo c" ) );
    JPanel charts = new JPanel( new GridLayout( 2, 1 ) );
    charts.add( cd ); charts.add( cr );
    JLabel info = new JLabel( blockSummary() );
    info.setBorder( BorderFactory.createEmptyBorder( 4, 10, 4, 10 ) );
    JPanel p = new JPanel( new BorderLayout() );
    p.add( charts, BorderLayout.CENTER ); p.add( info, BorderLayout.SOUTH );
    return p;
  }

  private String blockSummary() {
    csDifusividade.Result r = myR;
    double mean = 0, sd = 0; int m = 0;
    for( double d : r.dBlk ) if( !Double.isNaN( d ) ) { mean += d; m++; }
    mean = m > 0 ? mean / m : Double.NaN;
    for( double d : r.dBlk ) if( !Double.isNaN( d ) ) sd += ( d - mean ) * ( d - mean );
    sd = m > 1 ? Math.sqrt( sd / ( m - 1 ) ) : 0;
    double to = r.dTotal >= 0 ? r.azimuth : ( r.azimuth + 180 ) % 360, from = ( to + 180 ) % 360;
    StringBuilder sb = new StringBuilder( "<html>" );
    sb.append( String.format( Locale.US, "D = energia do f–k em %.3g–%.3g Hz e %.0f–%.0f m/s andando para +s (azimute %.0f°) menos a que anda para −s, dividida pela soma. "
        + "D = 0: simétrico; ±1: tudo de um lado. Faixa sombreada: |D| &lt; 0,3.<br>", r.p.f1, r.p.f2, r.p.cmin, r.p.cmax, r.azimuth ) );
    sb.append( String.format( Locale.US, "<b>Registro todo: D = %+.2f</b>; por bloco: média %+.2f, desvio %.2f (%d blocos). ", r.dTotal, mean, sd, m ) );
    if( Math.abs( r.dTotal ) < 0.3 ) sb.append( "<font color='#2f7d3c'>Campo aproximadamente simétrico ao longo da linha.</font>" );
    else sb.append( String.format( Locale.US, "<font color='#9b2f2f'>Mais energia andando para %.0f° (%s), isto é, chegando de ~%.0f° (%s).</font>", to, rumo( to ), from, rumo( from ) ) );
    if( m > 2 ) sb.append( sd < 0.1 ? " Estável ao longo do tempo." : " Varia ao longo do tempo: compare com o tempo (vento, ondas) e com atividades na área." );
    return sb.append( "</html>" ).toString();
  }

  //====================================================================
  // 4. Neighbour delays (clock / time alignment)
  //====================================================================
  private JComponent delayTab() {
    csDifusividade.Result r = myR;
    int m = r.excessT.length, nb = r.nblk;
    double[] thrA = new double[1];
    boolean[] sus = csDifusividade.suspects( r, thrA );
    double thr = thrA[0];
    String[] cols = new String[4 + nb];
    cols[0] = "Par"; cols[1] = "Δs (m)"; cols[2] = "atraso (ms)"; cols[3] = "excesso (ms)";
    for( int b = 0; b < nb; b++ ) cols[4 + b] = String.format( Locale.US, "bloco %d", b + 1 );
    Object[][] data = new Object[m][cols.length];
    for( int k = 0; k < m; k++ ) {
      data[k][0] = r.labels[r.order[k]] + " – " + r.labels[r.order[k + 1]] + ( sus[k] ? "  ⚠" : "" );
      data[k][1] = String.format( Locale.US, "%.0f", r.ds[k] );
      data[k][2] = ms( r.tauT[k] );
      data[k][3] = ms( r.excessT[k] );
      for( int b = 0; b < nb; b++ ) data[k][4 + b] = ms( r.excessB[b][k] );
    }
    DefaultTableModel model = new DefaultTableModel( data, cols ) { @Override public boolean isCellEditable( int a, int b ) { return false; } };
    JTable tab = new JTable( model );
    tab.setAutoResizeMode( nb > 8 ? JTable.AUTO_RESIZE_OFF : JTable.AUTO_RESIZE_ALL_COLUMNS );
    tab.getColumnModel().getColumn( 0 ).setPreferredWidth( 140 );
    DefaultTableCellRenderer rend = new DefaultTableCellRenderer() {
      @Override
      public Component getTableCellRendererComponent( JTable t, Object v, boolean sel, boolean foc, int row, int col ) {
        Component c = super.getTableCellRendererComponent( t, v, sel, foc, row, col );
        setHorizontalAlignment( col == 0 ? LEFT : RIGHT );
        Color bg = Color.white;
        if( col >= 3 ) {
          double e = col == 3 ? r.excessT[row] : r.excessB[col - 4][row];
          if( !Double.isNaN( e ) ) {
            double a = Math.min( 1.0, Math.abs( e ) / ( 2 * thr ) );
            bg = e > 0 ? new Color( 255, (int)( 255 - 120 * a ), (int)( 255 - 120 * a ) ) : new Color( (int)( 255 - 120 * a ), (int)( 255 - 90 * a ), 255 );
          }
        }
        if( !sel ) c.setBackground( bg );
        c.setFont( c.getFont().deriveFont( sus[row] ? Font.BOLD : Font.PLAIN ) );
        return c;
      }
    };
    for( int c = 0; c < cols.length; c++ ) tab.getColumnModel().getColumn( c ).setCellRenderer( rend );

    Chart ch = new Chart( "par de vizinhos (ordem ao longo da linha)", "excesso de atraso (ms)", myStatus );
    double[] x = new double[m], y = new double[m];
    double ym = 2 * thr * 1000;
    for( int k = 0; k < m; k++ ) { x[k] = k + 1; y[k] = r.excessT[k] * 1000; if( !Double.isNaN( y[k] ) ) ym = Math.max( ym, Math.abs( y[k] ) * 1.1 ); }
    ch.setRange( 0.5, m + 0.5, -ym, ym );
    ch.bandLo = -thr * 1000; ch.bandHi = thr * 1000;
    ch.hLines.add( 0.0 );
    for( int b = 0; b < nb; b++ ) {
      double[] yb = new double[m];
      for( int k = 0; k < m; k++ ) yb[k] = r.excessB[b][k] * 1000;
      ch.series.add( new Ser( x, yb, GREY, Ser.DOT_SMALL, b == 0 ? "cada bloco" : null ) );
    }
    ch.series.add( new Ser( x, y, Color.black, Ser.LINE_DOTS, "registro todo" ) );
    ch.labeler = xv -> {
      int k = (int)Math.round( xv ) - 1;
      return k >= 0 && k < m ? r.labels[r.order[k]] + "–" + r.labels[r.order[k + 1]] : "";
    };

    StringBuilder sb = new StringBuilder( "<html>" );
    sb.append( String.format( Locale.US, "Atraso de cada receptor em relação ao vizinho anterior, pela inclinação da fase da coerência em %.3g–%.3g Hz. "
        + "Excesso = atraso − (atraso de propagação mediano × Δs), com propagação aparente de %.0f m/s. "
        + "<b>Um erro de relógio ou de alinhamento dos tempos é constante</b>: aparece com o mesmo valor em todos os blocos (e em qualquer banda); "
        + "o campo de ondas muda de bloco para bloco.<br>", r.p.fd1, r.p.fd2, Math.abs( 1 / r.slownessT ) ) );
    StringBuilder list = new StringBuilder();
    for( int k = 0; k < m; k++ ) if( sus[k] ) {
      list.append( String.format( Locale.US, "%s–%s (%+.0f ms)  ", r.labels[r.order[k]], r.labels[r.order[k + 1]], r.excessT[k] * 1000 ) );
    }
    if( list.length() > 0 ) {
      sb.append( String.format( Locale.US, "<font color='#9b2f2f'><b>Suspeitos</b> (|excesso| &gt; %.0f ms no registro todo e mesmo sinal em ≥ 75%% dos blocos): %s</font><br>", thr * 1000, list ) );
      sb.append( "Um salto num par e o oposto no par seguinte = um receptor deslocado; um salto isolado = todos os receptores de um lado deslocados. Confira o tempo inicial dos arquivos desses receptores." );
    }
    else {
      sb.append( String.format( Locale.US, "<font color='#2f7d3c'>Nenhum par com excesso estável acima de %.0f ms.</font>", thr * 1000 ) );
    }
    sb.append( "</html>" );
    JLabel info = new JLabel( sb.toString() );
    info.setBorder( BorderFactory.createEmptyBorder( 4, 10, 4, 10 ) );
    JPanel p = new JPanel( new BorderLayout() );
    JSplitPane sp = new JSplitPane( JSplitPane.VERTICAL_SPLIT, ch, new JScrollPane( tab ) );
    sp.setResizeWeight( 0.5 );
    sp.setDividerLocation( 300 );
    p.add( sp, BorderLayout.CENTER ); p.add( info, BorderLayout.SOUTH );
    return p;
  }

  static String ms( double s ) { return Double.isNaN( s ) ? "" : String.format( Locale.US, "%+.0f", s * 1000 ); }

  //====================================================================
  private void exportCsv() {
    JFileChooser fc = new JFileChooser();
    fc.setDialogTitle( "Nome base dos arquivos CSV" );
    fc.setSelectedFile( new File( "difusividade" ) );
    if( fc.showSaveDialog( this ) != JFileChooser.APPROVE_OPTION ) return;
    String base = fc.getSelectedFile().getPath().replaceAll( "\\.csv$", "" );
    csDifusividade.Result r = myR;
    try {
      try( PrintWriter w = new PrintWriter( base + "_spac.csv", "UTF-8" ) ) {
        w.println( "f_Hz,c_ajuste_m_s,rms_J0,rms_cos,c_primeiro_zero_m_s,a,b,claro" );
        for( int g = 0; g < r.freqs.length; g++ ) w.println( String.format( Locale.US, "%.4f,%.1f,%.4f,%.4f,%.1f,%.3f,%.3f,%d",
            r.freqs[g], r.cFit[g], r.rmsJ0[g], r.rmsCos[g], r.cZero[g], r.ampA[g], r.ampB[g], clear( g ) ? 1 : 0 ) );
      }
      try( PrintWriter w = new PrintWriter( base + "_simetria.csv", "UTF-8" ) ) {
        w.println( "bloco,inicio_s,janelas,D,rms_J0,rms_cos" );
        for( int b = 0; b < r.nblk; b++ ) w.println( String.format( Locale.US, "%d,%.1f,%d,%.4f,%.4f,%.4f", b + 1, r.blkStart[b], r.blkWins[b], r.dBlk[b], r.rmsJ0Blk[b], r.rmsCosBlk[b] ) );
      }
      try( PrintWriter w = new PrintWriter( base + "_atrasos.csv", "UTF-8" ) ) {
        StringBuilder h = new StringBuilder( "receptor_a,receptor_b,ds_m,atraso_ms,excesso_ms" );
        for( int b = 0; b < r.nblk; b++ ) h.append( ",excesso_bloco" ).append( b + 1 ).append( "_ms" );
        w.println( h );
        for( int k = 0; k < r.excessT.length; k++ ) {
          StringBuilder l = new StringBuilder( String.format( Locale.US, "%s,%s,%.1f,%.1f,%.1f", r.labels[r.order[k]], r.labels[r.order[k + 1]], r.ds[k], r.tauT[k] * 1000, r.excessT[k] * 1000 ) );
          for( int b = 0; b < r.nblk; b++ ) l.append( String.format( Locale.US, ",%.1f", r.excessB[b][k] * 1000 ) );
          w.println( l );
        }
      }
      try( PrintWriter w = new PrintWriter( base + "_fk.csv", "UTF-8" ) ) {
        StringBuilder h = new StringBuilder( "k_1_km" );
        for( double f : r.freqs ) h.append( String.format( Locale.US, ",%.4f", f ) );
        w.println( h );
        for( int k = 0; k < r.ks.length; k++ ) {
          StringBuilder l = new StringBuilder( String.format( Locale.US, "%.5f", r.ks[k] * 1000 ) );
          for( int g = 0; g < r.freqs.length; g++ ) l.append( String.format( Locale.US, ",%.4f", r.fk[k][g] ) );
          w.println( l );
        }
      }
      myStatus.setText( "Salvos: " + new File( base ).getName() + "_spac/_simetria/_atrasos/_fk.csv" );
    }
    catch( IOException ex ) {
      JOptionPane.showMessageDialog( this, "Erro ao salvar:\n" + ex.getMessage() );
    }
  }

  private void savePng() {
    JFileChooser fc = new JFileChooser();
    fc.setSelectedFile( new File( "difusividade_aba" + ( myTabs.getSelectedIndex() + 1 ) + ".png" ) );
    if( fc.showSaveDialog( this ) != JFileChooser.APPROVE_OPTION ) return;
    File f = fc.getSelectedFile();
    if( !f.getName().toLowerCase().endsWith( ".png" ) ) f = new File( f.getPath() + ".png" );
    try {
      javax.imageio.ImageIO.write( render( (JComponent)myTabs.getSelectedComponent() ), "png", f );
      myStatus.setText( "Salvo: " + f.getName() );
    }
    catch( IOException ex ) {
      JOptionPane.showMessageDialog( this, "Erro ao salvar:\n" + ex.getMessage() );
    }
  }
  static BufferedImage render( JComponent c ) {
    BufferedImage bi = new BufferedImage( Math.max( 1, c.getWidth() ), Math.max( 1, c.getHeight() ), BufferedImage.TYPE_INT_RGB );
    Graphics2D g = bi.createGraphics();
    c.paint( g );
    g.dispose();
    return bi;
  }
  /** For tests: renders every tab to base_N.png */
  public static void snap( csDifusividade.Result r, String title, String base ) throws Exception {
    csDifusividadeFrame f = new csDifusividadeFrame( r, title );
    f.setVisible( true );
    for( int i = 0; i < f.myTabs.getTabCount(); i++ ) {
      final int ii = i;
      javax.swing.SwingUtilities.invokeAndWait( () -> f.myTabs.setSelectedIndex( ii ) );
      Thread.sleep( 700 );
      final BufferedImage[] im = new BufferedImage[1];
      javax.swing.SwingUtilities.invokeAndWait( () -> im[0] = render( (JComponent)f.getContentPane() ) );
      javax.imageio.ImageIO.write( im[0], "png", new File( base + "_" + ( i + 1 ) + ".png" ) );
    }
    f.dispose();
  }

  //====================================================================
  // Plot primitives
  //====================================================================
  static final class Ser {
    static final int LINE = 0, DASH = 1, DOT = 2, OPEN = 3, CROSS = 4, DOT_SMALL = 5, LINE_DOTS = 6;
    final double[] x, y; final Color c; final int kind; final String name;
    Ser( double[] x, double[] y, Color c, int kind, String name ) { this.x = x; this.y = y; this.c = c; this.kind = kind; this.name = name; }
    static Ser of( List<double[]> pts, Color c, int kind, String name ) {
      double[] x = new double[pts.size()], y = new double[pts.size()];
      for( int i = 0; i < x.length; i++ ) { x[i] = pts.get( i )[0]; y[i] = pts.get( i )[1]; }
      return new Ser( x, y, c, kind, name );
    }
  }

  interface XYMap { int xp( double x ); int yp( double y ); }

  static void drawSeries( Graphics2D g, List<Ser> list, XYMap m ) {
    Stroke s0 = g.getStroke();
    for( Ser s : list ) {
      g.setColor( s.c );
      if( s.kind == Ser.LINE || s.kind == Ser.DASH || s.kind == Ser.LINE_DOTS ) {
        g.setStroke( s.kind == Ser.DASH ? new BasicStroke( 1.4f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND, 1f, new float[]{ 6f, 4f }, 0f ) : new BasicStroke( 1.8f ) );
        int px = 0, py = 0; boolean has = false;
        for( int i = 0; i < s.x.length; i++ ) {
          if( Double.isNaN( s.y[i] ) || Double.isNaN( s.x[i] ) ) { has = false; continue; }
          int xx = m.xp( s.x[i] ), yy = m.yp( s.y[i] );
          if( has ) g.drawLine( px, py, xx, yy );
          px = xx; py = yy; has = true;
        }
        g.setStroke( s0 );
      }
      if( s.kind == Ser.LINE_DOTS || s.kind == Ser.DOT || s.kind == Ser.OPEN || s.kind == Ser.CROSS || s.kind == Ser.DOT_SMALL ) {
        int r = s.kind == Ser.DOT_SMALL ? 2 : 4;
        for( int i = 0; i < s.x.length; i++ ) {
          if( Double.isNaN( s.y[i] ) || Double.isNaN( s.x[i] ) ) continue;
          int xx = m.xp( s.x[i] ), yy = m.yp( s.y[i] );
          if( s.kind == Ser.OPEN ) { g.setStroke( new BasicStroke( 1.4f ) ); g.drawOval( xx - r, yy - r, 2 * r, 2 * r ); g.setStroke( s0 ); }
          else if( s.kind == Ser.CROSS ) { g.setStroke( new BasicStroke( 1.6f ) ); g.drawLine( xx - r, yy - r, xx + r, yy + r ); g.drawLine( xx - r, yy + r, xx + r, yy - r ); g.setStroke( s0 ); }
          else g.fillOval( xx - r, yy - r, 2 * r, 2 * r );
        }
      }
    }
  }

  static void drawLegend( Graphics2D g, List<Ser> list, int xr, int yt, boolean dark ) {
    FontMetrics fm = g.getFontMetrics();
    int w = 0, n = 0;
    for( Ser s : list ) if( s.name != null ) { w = Math.max( w, fm.stringWidth( s.name ) ); n++; }
    if( n == 0 ) return;
    int bw = w + 36, bh = n * 15 + 6, x0 = xr - bw - 6, y0 = yt + 6;
    g.setColor( dark ? new Color( 0, 0, 0, 150 ) : new Color( 255, 255, 255, 220 ) );
    g.fillRect( x0, y0, bw, bh );
    g.setColor( Color.gray ); g.drawRect( x0, y0, bw, bh );
    int y = y0 + 14;
    for( Ser s : list ) {
      if( s.name == null ) continue;
      List<Ser> one = new ArrayList<Ser>();
      one.add( new Ser( s.kind == Ser.LINE || s.kind == Ser.DASH || s.kind == Ser.LINE_DOTS ? new double[]{ 0, 1 } : new double[]{ 0.5 },
                        s.kind == Ser.LINE || s.kind == Ser.DASH || s.kind == Ser.LINE_DOTS ? new double[]{ 0, 0 } : new double[]{ 0 }, s.c, s.kind, null ) );
      final int yy = y - 4, xs = x0 + 6;
      drawSeries( g, one, new XYMap() { public int xp( double x ) { return xs + (int)Math.round( x * 20 ); } public int yp( double v ) { return yy; } } );
      g.setColor( dark ? Color.white : Color.darkGray );
      g.drawString( s.name, x0 + 30, y );
      y += 15;
    }
  }

  /** Simple x-y chart */
  static final class Chart extends JPanel implements XYMap {
    static final int ML = 66, MR = 18, MT = 12, MB = 42;
    double x0 = 0, x1 = 1, y0 = 0, y1 = 1, bandLo = Double.NaN, bandHi = Double.NaN;
    final List<Ser> series = new ArrayList<Ser>();
    final List<Double> hLines = new ArrayList<Double>();
    final String xl, yl;
    java.util.function.DoubleFunction<String> labeler;
    Chart( String xl, String yl, JLabel readout ) {
      this.xl = xl; this.yl = yl;
      setBackground( Color.white );
      setPreferredSize( new Dimension( 900, 300 ) );
      addMouseMotionListener( new MouseAdapter() {
        @Override public void mouseMoved( MouseEvent e ) {
          double x = x0 + ( e.getX() - ML ) / (double)Math.max( 1, getWidth() - ML - MR ) * ( x1 - x0 );
          double y = y1 - ( e.getY() - MT ) / (double)Math.max( 1, getHeight() - MT - MB ) * ( y1 - y0 );
          String lab = labeler != null ? labeler.apply( x ) + ": " : String.format( Locale.US, "x = %.4g, ", x );
          readout.setText( lab + String.format( Locale.US, "y = %.3g", y ) );
        }
      } );
    }
    void setRange( double a, double b, double c, double d ) { x0 = a; x1 = b; y0 = c; y1 = d; }
    public int xp( double x ) { return ML + (int)Math.round( ( x - x0 ) / ( x1 - x0 ) * ( getWidth() - ML - MR ) ); }
    public int yp( double y ) { return MT + (int)Math.round( ( y1 - y ) / ( y1 - y0 ) * ( getHeight() - MT - MB ) ); }
    @Override
    protected void paintComponent( Graphics g0 ) {
      super.paintComponent( g0 );
      Graphics2D g = (Graphics2D)g0;
      g.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
      g.setFont( g.getFont().deriveFont( Font.PLAIN, 11f ) );
      FontMetrics fm = g.getFontMetrics();
      int w = getWidth() - ML - MR, h = getHeight() - MT - MB;
      g.setColor( new Color( 250, 250, 250 ) ); g.fillRect( ML, MT, w, h );
      if( !Double.isNaN( bandLo ) ) { g.setColor( new Color( 225, 235, 240 ) ); int a = yp( bandHi ), b = yp( bandLo ); g.fillRect( ML, a, w, b - a ); }
      g.setColor( Color.gray );
      for( double t : csDispersionFC.ticks( y0, y1, Math.max( 2, h / 40 ) ) ) {
        int y = yp( t ); g.drawLine( ML - 4, y, ML, y );
        String s = csDispersionFC.fmt( t ); g.drawString( s, ML - 8 - fm.stringWidth( s ), y + 4 );
      }
      for( double t : csDispersionFC.ticks( x0, x1, Math.max( 2, w / 80 ) ) ) {
        int x = xp( t ); g.drawLine( x, MT + h, x, MT + h + 4 );
        String s = csDispersionFC.fmt( t ); g.drawString( s, x - fm.stringWidth( s ) / 2, MT + h + 16 );
      }
      g.setColor( Color.darkGray );
      for( double v : hLines ) g.drawLine( ML, yp( v ), ML + w, yp( v ) );
      g.drawString( xl, ML + w / 2 - fm.stringWidth( xl ) / 2, MT + h + 33 );
      Graphics2D gr = (Graphics2D)g.create(); gr.rotate( -Math.PI / 2 );
      gr.drawString( yl, -( MT + h / 2 ) - fm.stringWidth( yl ) / 2, 15 ); gr.dispose();
      Graphics2D gc = (Graphics2D)g.create(); gc.clipRect( ML, MT, w + 1, h + 1 );
      drawSeries( gc, series, this ); gc.dispose();
      g.setColor( Color.darkGray ); g.drawRect( ML, MT, w, h );
      drawLegend( g, series, ML + w, MT, false );
    }
  }

  /** Image v[ix][iy] on monotonic axes xs, ys (nearest sample), with overlays in data coordinates */
  static final class Heat extends JPanel implements XYMap {
    static final int ML = 70, MR = 70, MT = 12, MB = 42;
    final double[] xs, ys; final double[][] v; final double vmin, vmax;
    final String xl, yl;
    final List<Ser> series = new ArrayList<Ser>();
    final List<Double> vLines = new ArrayList<Double>();
    final int[] lut = csDispersionFC.colormapLut( "viridis" );
    java.util.function.BiFunction<Double, Double, String> reader;
    java.util.function.BiConsumer<Double, Double> click;
    private BufferedImage myImg; private int myW = -1, myH = -1;
    Heat( double[] xs, double[] ys, double[][] v, double vmin, double vmax, String xl, String yl, JLabel readout ) {
      this.xs = xs; this.ys = ys; this.v = v; this.vmin = vmin; this.vmax = vmax; this.xl = xl; this.yl = yl;
      setBackground( Color.white );
      setPreferredSize( new Dimension( 900, 360 ) );
      MouseAdapter ma = new MouseAdapter() {
        @Override public void mouseMoved( MouseEvent e ) {
          if( reader != null && inPlot( e ) ) readout.setText( reader.apply( xd( e.getX() ), yd( e.getY() ) ) );
        }
        @Override public void mouseClicked( MouseEvent e ) {
          if( click != null && inPlot( e ) ) click.accept( xd( e.getX() ), yd( e.getY() ) );
        }
      };
      addMouseListener( ma ); addMouseMotionListener( ma );
    }
    boolean inPlot( MouseEvent e ) { return e.getX() >= ML && e.getX() <= getWidth() - MR && e.getY() >= MT && e.getY() <= getHeight() - MB; }
    double xa() { return xs[0]; } double xb() { return xs[xs.length - 1]; }
    double ya() { return ys[0]; } double yb() { return ys[ys.length - 1]; }
    double xd( int px ) { return xa() + ( px - ML ) / (double)Math.max( 1, getWidth() - ML - MR ) * ( xb() - xa() ); }
    double yd( int py ) { return yb() - ( py - MT ) / (double)Math.max( 1, getHeight() - MT - MB ) * ( yb() - ya() ); }
    public int xp( double x ) { return ML + (int)Math.round( ( x - xa() ) / ( xb() - xa() ) * ( getWidth() - ML - MR ) ); }
    public int yp( double y ) { return MT + (int)Math.round( ( yb() - y ) / ( yb() - ya() ) * ( getHeight() - MT - MB ) ); }
    static int idx( double[] a, double x ) {
      int lo = 0, hi = a.length - 1;
      while( hi - lo > 1 ) { int m = ( lo + hi ) / 2; if( a[m] <= x ) lo = m; else hi = m; }
      return ( x - a[lo] <= a[hi] - x ) ? lo : hi;
    }
    @Override
    protected void paintComponent( Graphics g0 ) {
      super.paintComponent( g0 );
      Graphics2D g = (Graphics2D)g0;
      g.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
      g.setFont( g.getFont().deriveFont( Font.PLAIN, 11f ) );
      FontMetrics fm = g.getFontMetrics();
      int w = Math.max( 1, getWidth() - ML - MR ), h = Math.max( 1, getHeight() - MT - MB );
      if( myImg == null || w != myW || h != myH ) {
        myImg = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        int[] ix = new int[w];
        for( int px = 0; px < w; px++ ) ix[px] = idx( xs, xa() + ( px + 0.5 ) / w * ( xb() - xa() ) );
        for( int py = 0; py < h; py++ ) {
          int iy = idx( ys, yb() - ( py + 0.5 ) / h * ( yb() - ya() ) );
          for( int px = 0; px < w; px++ ) {
            double val = v[ix[px]][iy];
            int rgb = Double.isNaN( val ) ? 0x404040 : lut[(int)Math.max( 0, Math.min( 255, Math.round( 255 * ( val - vmin ) / ( vmax - vmin ) ) ) )];
            myImg.setRGB( px, py, rgb );
          }
        }
        myW = w; myH = h;
      }
      g.drawImage( myImg, ML, MT, null );
      g.setColor( Color.gray );
      for( double t : csDispersionFC.ticks( ya(), yb(), Math.max( 2, h / 40 ) ) ) {
        int y = yp( t ); g.drawLine( ML - 4, y, ML, y );
        String s = csDispersionFC.fmt( t ); g.drawString( s, ML - 8 - fm.stringWidth( s ), y + 4 );
      }
      for( double t : csDispersionFC.ticks( xa(), xb(), Math.max( 2, w / 80 ) ) ) {
        int x = xp( t ); g.drawLine( x, MT + h, x, MT + h + 4 );
        String s = csDispersionFC.fmt( t ); g.drawString( s, x - fm.stringWidth( s ) / 2, MT + h + 16 );
      }
      g.setColor( Color.darkGray );
      g.drawString( xl, ML + w / 2 - fm.stringWidth( xl ) / 2, MT + h + 33 );
      Graphics2D gr = (Graphics2D)g.create(); gr.rotate( -Math.PI / 2 );
      gr.drawString( yl, -( MT + h / 2 ) - fm.stringWidth( yl ) / 2, 15 ); gr.dispose();
      Graphics2D gc = (Graphics2D)g.create(); gc.clipRect( ML, MT, w + 1, h + 1 );
      gc.setColor( new Color( 255, 255, 255, 150 ) );
      for( double f : vLines ) gc.drawLine( xp( f ), MT, xp( f ), MT + h );
      drawSeries( gc, series, this ); gc.dispose();
      g.setColor( Color.darkGray ); g.drawRect( ML, MT, w, h );
      drawLegend( g, series, ML + w, MT, true );
      // colour bar
      int bx = ML + w + 14, bw = 12;
      for( int py = 0; py < h; py++ ) { g.setColor( new Color( lut[(int)Math.round( 255.0 * ( h - 1 - py ) / Math.max( 1, h - 1 ) )] ) ); g.drawLine( bx, MT + py, bx + bw, MT + py ); }
      g.setColor( Color.darkGray ); g.drawRect( bx, MT, bw, h );
      g.drawString( csDispersionFC.fmt( vmax ), bx + bw + 3, MT + 10 ); g.drawString( csDispersionFC.fmt( vmin ), bx + bw + 3, MT + h );
    }
  }
}
