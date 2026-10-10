/* f-c dispersion image (phase-shift method) of a gather of correlations (VSG)
   M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import javax.imageio.ImageIO;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.filechooser.FileNameExtensionFilter;

/**
 * Phase-velocity vs frequency (f-c) image of a gather of correlations, e.g. a virtual shot gather (VSG):
 * <pre>
 *   o_x(t)  = one side of the correlation at offset x (lag 0 ... tmax), or both sides summed
 *   U(f,x)  = O(f,x) / |O(f,x)|                (spectrum with unit amplitude: only the phase is used)
 *   E(f,c)  = | sum_x U(f,x) exp( i 2 pi f x / c ) | / N        (0..1; 1 = all traces in phase)
 * </pre>
 * (Park et al., 1998, "Imaging dispersion curves of surface waves on multi-channel record").
 * A wave arriving at t = x/c(f) stacks in phase at that velocity, independent of the amplitudes.
 */
public final class csDispersionFC {
  public static final int SIDE_CAUSAL = 0, SIDE_ACAUSAL = 1, SIDE_SUM = 2;
  public static final String[] SIDES = { "lags positivos (causal)", "lags negativos (acausal)", "soma dos dois lados" };
  public static final int WIN_ALL = 0, WIN_FAST = 1, WIN_SLOW = 2;
  public static final String[] WINDOWS = { "traço todo", "rápida: t < x/v + folga", "lenta: t > x/v + folga", "rápida e lenta (duas imagens)" };

  /** Weight of the time window at lag t [s] for offset x [m] (smooth tanh boundary of half-width p.edge) */
  static double windowWeight( Params p, double t, double x ) {
    if( p.timeWindow == WIN_ALL ) return 1.0;
    double z = ( t - ( x / p.vcut + p.pad ) ) / Math.max( 1.0e-3, p.edge );
    double slow = 0.5 + 0.5 * Math.tanh( z );
    return ( p.timeWindow == WIN_SLOW ) ? slow : 1.0 - slow;
  }

  /** Parameters (s, Hz, m, m/s) */
  public static final class Params {
    public int    lag0 = 0;            // sample of lag 0 in the input traces
    public int    side = SIDE_SUM;     // which side of the correlation (ignored if the input has only lags >= 0)
    public boolean oneSided = false;   // input already has only lags 0..max (e.g. VSG with both sides summed)
    public double tmax = 0.0;          // maximum lag used [s], 0 = all
    public double fmin = 0.1, fmax = 2.5;
    public double cmin = 150, cmax = 2500;
    public int    nc = 300;
    public double dmin = 0, dmax = 1.0e30;
    public int    timeWindow = WIN_ALL;  // WIN_ALL, WIN_FAST (t < x/vcut + pad) or WIN_SLOW (t > x/vcut + pad)
    public double vcut = 1200, pad = 2.0, edge = 0.5;   // [m/s], [s], half-width of the smooth boundary [s]
  }

  /** Result */
  public static final class Image {
    public double[] freqs;     // [nf]
    public double[] vels;      // [nc]
    public double[][] E;       // [nc][nf], 0..1
    public int nUsed;          // traces used
    public double dmin, dmax;  // offset range actually used [m]
    public double dx;          // typical offset spacing [m] (spatial aliasing: c < 2 dx f)
    public int nLagSamples;
    public double dt;
  }

  /**
   * @param traces correlations [ntr][ns]
   * @param offs   offset of each trace [m] (absolute value is used; NaN: trace skipped)
   * @param dt     sample interval [s]
   */
  public static Image compute( float[][] traces, double[] offs, double dt, Params p ) {
    int ntr = traces.length;
    int ns = traces[0].length;
    int n;                                       // samples per one-sided trace
    if( p.oneSided ) n = ns - p.lag0;
    else if( p.side == SIDE_CAUSAL ) n = ns - p.lag0;
    else if( p.side == SIDE_ACAUSAL ) n = p.lag0 + 1;
    else n = Math.min( ns - p.lag0, p.lag0 + 1 );
    if( p.tmax > 0 ) n = Math.min( n, (int)Math.round( p.tmax / dt ) + 1 );
    if( n < 4 ) throw new IllegalArgumentException( "Poucas amostras depois do lag 0 (" + n + "): confira o tempo do lag 0" );
    int nfft = 1;
    while( nfft < 2 * n ) nfft *= 2;
    double df = 1.0 / ( nfft * dt );
    double fNyq = 0.5 / dt;
    double fmax = Math.min( p.fmax > 0 ? p.fmax : fNyq, fNyq );
    int k0 = Math.max( 1, (int)Math.ceil( p.fmin / df ) ), k1 = Math.min( nfft / 2, (int)Math.floor( fmax / df ) );
    if( k1 < k0 ) throw new IllegalArgumentException( "Banda de frequência vazia" );
    int nf = k1 - k0 + 1;
    int nc = Math.max( 2, p.nc );

    Image img = new Image();
    img.dt = dt;
    img.nLagSamples = n;
    img.freqs = new double[nf];
    for( int j = 0; j < nf; j++ ) img.freqs[j] = ( k0 + j ) * df;
    img.vels = new double[nc];
    for( int i = 0; i < nc; i++ ) img.vels[i] = p.cmin + i * ( p.cmax - p.cmin ) / ( nc - 1 );
    double[][] accRe = new double[nc][nf], accIm = new double[nc][nf];

    double[] taper = new double[n];               // cosine taper over the last 10% of the lags
    int ne = Math.max( 1, n / 10 );
    for( int j = 0; j < n; j++ ) taper[j] = ( j < n - ne ) ? 1.0 : 0.5 + 0.5 * Math.cos( Math.PI * ( j - ( n - ne ) ) / ne );

    double[] re = new double[nfft], im = new double[nfft];
    double[] ur = new double[nf], ui = new double[nf];
    List<Double> used = new ArrayList<Double>();
    for( int it = 0; it < ntr; it++ ) {
      double x = Math.abs( offs[it] );
      if( Double.isNaN( x ) || x < p.dmin || x > p.dmax ) continue;
      float[] s = traces[it];
      Arrays.fill( re, 0.0 );
      Arrays.fill( im, 0.0 );
      boolean any = false;
      for( int j = 0; j < n; j++ ) {
        double v;
        if( p.oneSided || p.side == SIDE_CAUSAL ) v = s[p.lag0 + j];
        else if( p.side == SIDE_ACAUSAL ) v = s[p.lag0 - j];
        else v = s[p.lag0 + j] + s[p.lag0 - j];
        re[j] = v * taper[j] * windowWeight( p, j * dt, x );
        if( v != 0.0 ) any = true;
      }
      if( !any ) continue;                        // dead trace (e.g. no reference)
      csClockXcorr.fft( re, im );
      for( int j = 0; j < nf; j++ ) {
        double a = Math.hypot( re[k0 + j], im[k0 + j] );
        ur[j] = ( a > 0 ) ? re[k0 + j] / a : 0.0;
        ui[j] = ( a > 0 ) ? im[k0 + j] / a : 0.0;
      }
      // FFT convention exp(-i w t): a wave at t0 = x/c has phase -w t0 -> multiply by exp(+i w x / c)
      for( int ic = 0; ic < nc; ic++ ) {
        double tau = x / img.vels[ic];
        double[] ar = accRe[ic], ai = accIm[ic];
        for( int j = 0; j < nf; j++ ) {
          double ph = 2.0 * Math.PI * img.freqs[j] * tau;
          double c = Math.cos( ph ), sn = Math.sin( ph );
          ar[j] += ur[j] * c - ui[j] * sn;
          ai[j] += ur[j] * sn + ui[j] * c;
        }
      }
      used.add( x );
    }
    img.nUsed = used.size();
    if( img.nUsed < 2 ) throw new IllegalArgumentException( "Só " + img.nUsed + " traço(s) com offset na faixa pedida: são necessários pelo menos 2" );
    img.E = new double[nc][nf];
    for( int ic = 0; ic < nc; ic++ ) {
      for( int j = 0; j < nf; j++ ) img.E[ic][j] = Math.hypot( accRe[ic][j], accIm[ic][j] ) / img.nUsed;
    }
    double[] xs = new double[used.size()];
    for( int i = 0; i < xs.length; i++ ) xs[i] = used.get( i );
    Arrays.sort( xs );
    img.dmin = xs[0];
    img.dmax = xs[xs.length - 1];
    List<Double> dd = new ArrayList<Double>();
    for( int i = 1; i < xs.length; i++ ) if( xs[i] - xs[i - 1] > 1.0 ) dd.add( xs[i] - xs[i - 1] );
    if( dd.isEmpty() ) img.dx = 0;
    else {
      double[] d = new double[dd.size()];
      for( int i = 0; i < d.length; i++ ) d[i] = dd.get( i );
      Arrays.sort( d );
      img.dx = d[d.length / 2];
    }
    return img;
  }

  /** Index of the velocity with maximum E at each frequency */
  public static int[] maxPerFrequency( Image img ) {
    int nf = img.freqs.length;
    int[] r = new int[nf];
    for( int j = 0; j < nf; j++ ) {
      int best = 0;
      for( int ic = 1; ic < img.vels.length; ic++ ) if( img.E[ic][j] > img.E[best][j] ) best = ic;
      r[j] = best;
    }
    return r;
  }

  /** Parses "lag 0 em 40000 ms" from a VSG pane title; returns the lag 0 in ms, or -1 */
  public static double lag0FromTitle( String title ) {
    if( title == null ) return -1;
    java.util.regex.Matcher m = java.util.regex.Pattern.compile( "lag 0 em ([0-9.]+) ms" ).matcher( title );
    if( m.find() ) {
      try { return Double.parseDouble( m.group( 1 ) ); } catch( NumberFormatException e ) { return -1; }
    }
    return -1;
  }

  //====================================================================
  // Display
  //====================================================================
  static final String[] COLORMAPS = { "viridis", "turbo", "jet", "cinza invertido" };
  private static final double[][] CM_VIRIDIS = {
    {0,.267,.005,.329},{.1,.283,.141,.458},{.2,.254,.265,.530},{.3,.207,.372,.553},{.4,.164,.471,.558},{.5,.128,.567,.551},
    {.6,.135,.659,.518},{.7,.267,.749,.441},{.8,.478,.821,.318},{.9,.741,.873,.150},{1,.993,.906,.144} };
  private static final double[][] CM_JET = { {0,0,0,.5},{.125,0,0,1},{.375,0,1,1},{.625,1,1,0},{.875,1,0,0},{1,.5,0,0} };
  private static final double[][] CM_GRAY_INV = { {0,1,1,1},{1,0,0,0} };

  static int[] colormapLut( String name ) {
    int[] lut = new int[256];
    for( int i = 0; i < 256; i++ ) {
      double v = i / 255.0, r, g, b;
      if( "turbo".equals( name ) ) {
        r = 0.13572138 + v*(4.61539260 + v*(-42.66032258 + v*(132.13108234 + v*(-152.94239396 + v*59.28637943))));
        g = 0.09140261 + v*(2.19418839 + v*(4.84296658 + v*(-14.18503333 + v*(4.27729857 + v*2.82956604))));
        b = 0.10667330 + v*(12.64194608 + v*(-60.58204836 + v*(110.36276771 + v*(-89.90310912 + v*27.34824973))));
      }
      else {
        double[][] cm = "jet".equals( name ) ? CM_JET : "cinza invertido".equals( name ) ? CM_GRAY_INV : CM_VIRIDIS;
        int k = 0;
        while( k < cm.length - 2 && v > cm[k + 1][0] ) k++;
        double t = Math.max( 0, Math.min( 1, ( v - cm[k][0] ) / ( cm[k + 1][0] - cm[k][0] ) ) );
        r = cm[k][1] + t * ( cm[k + 1][1] - cm[k][1] );
        g = cm[k][2] + t * ( cm[k + 1][2] - cm[k][2] );
        b = cm[k][3] + t * ( cm[k + 1][3] - cm[k][3] );
      }
      int ri = (int)Math.round( 255 * Math.max( 0, Math.min( 1, r ) ) );
      int gi = (int)Math.round( 255 * Math.max( 0, Math.min( 1, g ) ) );
      int bi = (int)Math.round( 255 * Math.max( 0, Math.min( 1, b ) ) );
      lut[i] = ( ri << 16 ) | ( gi << 8 ) | bi;
    }
    return lut;
  }

  /** Window with the f-c image: cursor readout, maximum per frequency, aliasing line, picks, PNG/CSV export */
  public static final class Frame extends JFrame {
    private static final int ML = 70, MR = 80, MT = 14, MB = 46;
    private final Image myImg;
    private final String myTitle;
    private final int[] myMax;
    private final List<double[]> myPicks = new ArrayList<double[]>();   // {f, c}
    private final JLabel myCursor = new JLabel( " " );
    private final JCheckBox myNormF = new JCheckBox( "Normalizar por frequência", true );
    private final JCheckBox myShowMax = new JCheckBox( "Máximo por frequência", true );
    private final JCheckBox myShowAlias = new JCheckBox( "Limite de aliasing", true );
    private final JComboBox<String> myCmap = new JComboBox<String>( COLORMAPS );
    private int[] myLut = colormapLut( "viridis" );
    private final Plot myPlot = new Plot();

    public Frame( Image img, String title, String info ) {
      super( "Dispersão f-c — " + title );
      myImg = img;
      myTitle = title;
      myMax = maxPerFrequency( img );
      setDefaultCloseOperation( DISPOSE_ON_CLOSE );
      JLabel lInfo = new JLabel( "<html>" + info + "</html>" );
      lInfo.setBorder( BorderFactory.createEmptyBorder( 3, 8, 3, 8 ) );
      JPanel row1 = new JPanel( new FlowLayout( FlowLayout.LEFT, 6, 1 ) );
      row1.add( myNormF ); row1.add( myShowMax ); row1.add( myShowAlias );
      row1.add( new JLabel( "Cores:" ) ); row1.add( myCmap );
      JButton bPng = new JButton( "Salvar PNG..." ), bCsv = new JButton( "Exportar CSV..." ), bClr = new JButton( "Limpar picks" );
      JPanel row2 = new JPanel( new FlowLayout( FlowLayout.LEFT, 6, 1 ) );
      row2.add( bPng ); row2.add( bCsv ); row2.add( bClr );
      row2.add( new JLabel( "  Clique esquerdo: pick   direito: apaga pick  |" ) );
      row2.add( myCursor );
      JPanel top = new JPanel();
      top.setLayout( new BoxLayout( top, BoxLayout.Y_AXIS ) );
      for( javax.swing.JComponent c : new javax.swing.JComponent[]{ lInfo, row1, row2 } ) { c.setAlignmentX( 0f ); top.add( c ); }
      getContentPane().setLayout( new BorderLayout() );
      getContentPane().add( top, BorderLayout.NORTH );
      getContentPane().add( myPlot, BorderLayout.CENTER );
      myPlot.setPreferredSize( new Dimension( 900, 600 ) );
      myNormF.addActionListener( e -> { myPlot.myKey = ""; myPlot.repaint(); } );
      myShowMax.addActionListener( e -> myPlot.repaint() );
      myShowAlias.addActionListener( e -> myPlot.repaint() );
      myShowAlias.setEnabled( img.dx > 0 );
      myShowAlias.setToolTipText( String.format( Locale.US, "c = 2·Δx·f, Δx = %.0f m (espaçamento típico dos offsets). Abaixo da linha a velocidade é ambígua", img.dx ) );
      myCmap.addActionListener( e -> { myLut = colormapLut( (String)myCmap.getSelectedItem() ); myPlot.myKey = ""; myPlot.repaint(); } );
      bPng.addActionListener( e -> savePng() );
      bCsv.addActionListener( e -> saveCsv() );
      bClr.addActionListener( e -> { myPicks.clear(); myPlot.repaint(); } );
      pack();
    }

    double fA() { return myImg.freqs[0]; }
    double fB() { return myImg.freqs[myImg.freqs.length - 1]; }
    double cA() { return myImg.vels[0]; }
    double cB() { return myImg.vels[myImg.vels.length - 1]; }

    /** E at (f, c), bilinear; normalized per frequency if selected */
    double value( double f, double c ) {
      int nf = myImg.freqs.length, nc = myImg.vels.length;
      double xf = ( f - fA() ) / ( fB() - fA() ) * ( nf - 1 ), xc = ( c - cA() ) / ( cB() - cA() ) * ( nc - 1 );
      int j0 = Math.max( 0, Math.min( nf - 2, (int)Math.floor( xf ) ) ), i0 = Math.max( 0, Math.min( nc - 2, (int)Math.floor( xc ) ) );
      double wf = Math.max( 0, Math.min( 1, xf - j0 ) ), wc = Math.max( 0, Math.min( 1, xc - i0 ) );
      double v0 = norm( i0, j0 ) * ( 1 - wf ) + norm( i0, j0 + 1 ) * wf;
      double v1 = norm( i0 + 1, j0 ) * ( 1 - wf ) + norm( i0 + 1, j0 + 1 ) * wf;
      return v0 * ( 1 - wc ) + v1 * wc;
    }
    private double norm( int ic, int j ) {
      double v = myImg.E[ic][j];
      if( myNormF.isSelected() ) {
        double m = myImg.E[myMax[j]][j];
        return m > 0 ? v / m : 0;
      }
      return v;
    }
    private double globalMax() {
      double m = 0;
      for( double[] r : myImg.E ) for( double v : r ) m = Math.max( m, v );
      return m > 0 ? m : 1;
    }

    private void savePng() {
      JFileChooser fc = new JFileChooser();
      fc.setFileFilter( new FileNameExtensionFilter( "PNG (*.png)", "png" ) );
      fc.setSelectedFile( new File( "dispersao_fc.png" ) );
      if( fc.showSaveDialog( this ) != JFileChooser.APPROVE_OPTION ) return;
      File f = fc.getSelectedFile();
      if( !f.getName().toLowerCase().endsWith( ".png" ) ) f = new File( f.getPath() + ".png" );
      BufferedImage img = new BufferedImage( myPlot.getWidth(), myPlot.getHeight(), BufferedImage.TYPE_INT_RGB );
      Graphics2D g = img.createGraphics();
      myPlot.paint( g );
      g.dispose();
      try { ImageIO.write( img, "png", f ); myCursor.setText( "Salvo: " + f.getName() ); }
      catch( IOException ex ) { JOptionPane.showMessageDialog( this, "Erro ao salvar PNG:\n" + ex.getMessage() ); }
    }

    /** CSV: picks (if any), the maximum per frequency, and the whole matrix in a second file */
    private void saveCsv() {
      JFileChooser fc = new JFileChooser();
      fc.setFileFilter( new FileNameExtensionFilter( "CSV (*.csv)", "csv" ) );
      fc.setSelectedFile( new File( "dispersao_fc.csv" ) );
      if( fc.showSaveDialog( this ) != JFileChooser.APPROVE_OPTION ) return;
      File f = fc.getSelectedFile();
      if( !f.getName().toLowerCase().endsWith( ".csv" ) ) f = new File( f.getPath() + ".csv" );
      File fm = new File( f.getPath().replaceAll( "\\.csv$", "" ) + "_matriz.csv" );
      try( PrintWriter w = new PrintWriter( f, "UTF-8" ) ) {
        w.println( "# " + myTitle );
        w.println( "# tipo,frequencia_Hz,velocidade_fase_m_s,E,E_max_na_frequencia" );
        List<double[]> picks = new ArrayList<double[]>( myPicks );
        picks.sort( ( a, b ) -> Double.compare( a[0], b[0] ) );
        for( double[] pk : picks ) {
          int j = nearestF( pk[0] );
          w.println( String.format( Locale.US, "pick,%.5f,%.1f,%.4f,%.4f", pk[0], pk[1], rawValue( pk[0], pk[1] ), myImg.E[myMax[j]][j] ) );
        }
        for( int j = 0; j < myImg.freqs.length; j++ ) {
          w.println( String.format( Locale.US, "maximo,%.5f,%.1f,%.4f,%.4f", myImg.freqs[j], myImg.vels[myMax[j]], myImg.E[myMax[j]][j], myImg.E[myMax[j]][j] ) );
        }
      }
      catch( IOException ex ) { JOptionPane.showMessageDialog( this, "Erro ao salvar CSV:\n" + ex.getMessage() ); return; }
      try( PrintWriter w = new PrintWriter( fm, "UTF-8" ) ) {
        StringBuilder sb = new StringBuilder( "velocidade_m_s\\frequencia_Hz" );
        for( double fr : myImg.freqs ) sb.append( String.format( Locale.US, ",%.5f", fr ) );
        w.println( sb );
        for( int ic = 0; ic < myImg.vels.length; ic++ ) {
          sb = new StringBuilder( String.format( Locale.US, "%.1f", myImg.vels[ic] ) );
          for( int j = 0; j < myImg.freqs.length; j++ ) sb.append( String.format( Locale.US, ",%.5f", myImg.E[ic][j] ) );
          w.println( sb );
        }
      }
      catch( IOException ex ) { JOptionPane.showMessageDialog( this, "Erro ao salvar a matriz:\n" + ex.getMessage() ); return; }
      myCursor.setText( "Salvo: " + f.getName() + " e " + fm.getName() );
    }
    private int nearestF( double f ) {
      int j = (int)Math.round( ( f - fA() ) / ( fB() - fA() ) * ( myImg.freqs.length - 1 ) );
      return Math.max( 0, Math.min( myImg.freqs.length - 1, j ) );
    }
    private double rawValue( double f, double c ) {
      boolean n = myNormF.isSelected();
      myNormF.setSelected( false );
      double v = value( f, c );
      myNormF.setSelected( n );
      return v;
    }

    //------------------------------------------------------------------
    final class Plot extends JPanel {
      private BufferedImage myBuf;
      String myKey = "";
      Plot() {
        setBackground( Color.white );
        MouseAdapter ma = new MouseAdapter() {
          @Override public void mouseMoved( MouseEvent e ) { cursor( e.getX(), e.getY() ); }
          @Override public void mouseExited( MouseEvent e ) { myCursor.setText( " " ); }
          @Override public void mouseClicked( MouseEvent e ) { click( e ); }
        };
        addMouseListener( ma );
        addMouseMotionListener( ma );
      }
      int pw() { return Math.max( 20, getWidth() - ML - MR ); }
      int ph() { return Math.max( 20, getHeight() - MT - MB ); }
      double fAt( int x ) { return fA() + ( x - ML + 0.5 ) / pw() * ( fB() - fA() ); }
      double cAt( int y ) { return cB() - ( y - MT + 0.5 ) / ph() * ( cB() - cA() ); }
      int xOf( double f ) { return ML + (int)Math.round( ( f - fA() ) / ( fB() - fA() ) * pw() ); }
      int yOf( double c ) { return MT + (int)Math.round( ( cB() - c ) / ( cB() - cA() ) * ph() ); }
      boolean inPlot( int x, int y ) { return x >= ML && x < ML + pw() && y >= MT && y < MT + ph(); }

      private void cursor( int x, int y ) {
        if( !inPlot( x, y ) ) { myCursor.setText( " " ); return; }
        double f = fAt( x ), c = cAt( y );
        int j = nearestF( f );
        myCursor.setText( String.format( Locale.US, "f = %.3f Hz   c = %.0f m/s   E = %.3f   (máximo nesta f: %.0f m/s, E = %.3f;  λ = %.0f m)",
            f, c, rawValue( f, c ), myImg.vels[myMax[j]], myImg.E[myMax[j]][j], c / f ) );
      }
      private void click( MouseEvent e ) {
        if( !inPlot( e.getX(), e.getY() ) ) return;
        if( javax.swing.SwingUtilities.isRightMouseButton( e ) ) {
          int best = -1; double bd = 15 * 15;
          for( int i = 0; i < myPicks.size(); i++ ) {
            double dx = xOf( myPicks.get( i )[0] ) - e.getX(), dy = yOf( myPicks.get( i )[1] ) - e.getY();
            if( dx*dx + dy*dy < bd ) { bd = dx*dx + dy*dy; best = i; }
          }
          if( best >= 0 ) myPicks.remove( best );
        }
        else {
          myPicks.add( new double[]{ fAt( e.getX() ), cAt( e.getY() ) } );
        }
        repaint();
      }

      @Override
      protected void paintComponent( Graphics g0 ) {
        super.paintComponent( g0 );
        Graphics2D g = (Graphics2D)g0;
        g.setColor( Color.white );
        g.fillRect( 0, 0, getWidth(), getHeight() );
        g.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
        g.setFont( g.getFont().deriveFont( Font.PLAIN, 11f ) );
        FontMetrics fm = g.getFontMetrics();
        int pw = pw(), ph = ph();
        String key = pw + "x" + ph + myNormF.isSelected() + myCmap.getSelectedItem();
        if( myBuf == null || !key.equals( myKey ) ) {
          myBuf = new BufferedImage( pw, ph, BufferedImage.TYPE_INT_RGB );
          double scale = myNormF.isSelected() ? 1.0 : globalMax();
          for( int x = 0; x < pw; x++ ) {
            double f = fAt( ML + x );
            for( int y = 0; y < ph; y++ ) {
              double v = value( f, cAt( MT + y ) ) / scale;
              int k = (int)Math.round( 255 * Math.max( 0, Math.min( 1, v ) ) );
              myBuf.setRGB( x, y, myLut[k] );
            }
          }
          myKey = key;
        }
        g.drawImage( myBuf, ML, MT, null );
        java.awt.Shape clip = g.getClip();
        g.clipRect( ML, MT, pw, ph );
        if( myShowAlias.isSelected() && myImg.dx > 0 ) {
          g.setColor( new Color( 255, 60, 60 ) );
          g.setStroke( new BasicStroke( 1.5f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10f, new float[]{ 6f, 4f }, 0f ) );
          g.drawLine( xOf( fA() ), yOf( 2 * myImg.dx * fA() ), xOf( fB() ), yOf( 2 * myImg.dx * fB() ) );
        }
        if( myShowMax.isSelected() ) {
          g.setColor( Color.white );
          for( int j = 0; j < myImg.freqs.length; j++ ) {
            int x = xOf( myImg.freqs[j] ), y = yOf( myImg.vels[myMax[j]] );
            g.fillOval( x - 2, y - 2, 4, 4 );
          }
        }
        g.setStroke( new BasicStroke( 1.5f ) );
        for( double[] pk : myPicks ) {
          int x = xOf( pk[0] ), y = yOf( pk[1] );
          g.setColor( Color.black ); g.drawOval( x - 5, y - 5, 10, 10 );
          g.setColor( Color.magenta ); g.drawLine( x - 4, y, x + 4, y ); g.drawLine( x, y - 4, x, y + 4 );
        }
        g.setClip( clip );
        g.setStroke( new BasicStroke( 1f ) );
        g.setColor( Color.darkGray );
        g.drawRect( ML, MT, pw, ph );
        for( double f : ticks( fA(), fB(), Math.max( 2, pw / 80 ) ) ) {
          int x = xOf( f );
          g.drawLine( x, MT + ph, x, MT + ph + 4 );
          String s = fmt( f );
          g.drawString( s, x - fm.stringWidth( s ) / 2, MT + ph + 5 + fm.getAscent() );
        }
        String xl = "Frequência [Hz]";
        g.drawString( xl, ML + pw / 2 - fm.stringWidth( xl ) / 2, MT + ph + 8 + 2 * fm.getHeight() );
        for( double c : ticks( cA(), cB(), Math.max( 2, ph / 40 ) ) ) {
          int y = yOf( c );
          g.drawLine( ML - 4, y, ML, y );
          String s = fmt( c );
          g.drawString( s, ML - 6 - fm.stringWidth( s ), y + fm.getAscent() / 2 - 1 );
        }
        Graphics2D gr = (Graphics2D)g.create();
        gr.rotate( -Math.PI / 2 );
        String yl = "Velocidade de fase [m/s]";
        gr.drawString( yl, -( MT + ph / 2 ) - fm.stringWidth( yl ) / 2, 14 );
        gr.dispose();
        int cbx = ML + pw + 16, cbw = 14;
        for( int y = 0; y < ph; y++ ) {
          g.setColor( new Color( myLut[(int)Math.round( 255.0 * ( ph - 1 - y ) / Math.max( 1, ph - 1 ) )] ) );
          g.drawLine( cbx, MT + y, cbx + cbw, MT + y );
        }
        g.setColor( Color.darkGray );
        g.drawRect( cbx, MT, cbw, ph );
        double top = myNormF.isSelected() ? 1.0 : globalMax();
        for( int i = 0; i <= 4; i++ ) {
          int y = MT + ph - (int)Math.round( i / 4.0 * ph );
          g.drawLine( cbx + cbw, y, cbx + cbw + 3, y );
          g.drawString( String.format( Locale.US, "%.2f", top * i / 4.0 ), cbx + cbw + 5, y + fm.getAscent() / 2 - 1 );
        }
        g.drawString( "E", cbx + 3, MT - 2 );
      }
    }
  }

  static double[] ticks( double a, double b, int maxTicks ) {
    double raw = ( b - a ) / maxTicks;
    if( !( raw > 0 ) ) return new double[]{ a };
    double mag = Math.pow( 10, Math.floor( Math.log10( raw ) ) ), step = 10 * mag;
    for( double m : new double[]{ 1, 2, 2.5, 5, 10 } ) if( m * mag >= raw ) { step = m * mag; break; }
    List<Double> l = new ArrayList<Double>();
    for( double v = Math.ceil( a / step - 1e-9 ) * step; v <= b + 1e-9 * step; v += step ) l.add( v );
    double[] r = new double[l.size()];
    for( int i = 0; i < r.length; i++ ) r[i] = l.get( i );
    return r;
  }
  static String fmt( double v ) {
    if( Math.abs( v - Math.rint( v ) ) < 1e-9 ) return String.valueOf( (long)Math.rint( v ) );
    String s = String.format( Locale.US, "%.3f", v );
    return s.replaceAll( "0+$", "" ).replaceAll( "\\.$", "" );
  }
}
