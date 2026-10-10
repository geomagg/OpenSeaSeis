/* Causal / acausal asymmetry of noise cross-correlations (Interferometria)
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
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.filechooser.FileNameExtensionFilter;

/**
 * Asymmetry of two-sided correlations C(tau), tau = -nl..nl (lag 0 at sample nl):
 * <pre>
 *   A = (E+ - E-) / (E+ + E-),   E+- = sum over the arrival window of C(+-tau)^2
 *   arrival window: x/cmax - pad  &lt;=  tau  &lt;=  x/cmin + pad   (whole side if x is unknown)
 * </pre>
 * A &gt; 0: energy went from the virtual source to the trace; A &lt; 0: from the trace to the virtual source.
 * With the signed position s along the line (s &gt; 0 on one side of the source), D = A·sign(s) tells the direction
 * along the line in which the energy travels (D &gt; 0: towards +s). Diffuse field: A ~ 0.
 */
public final class csAssimetria {

  public static final class Params {
    public double cmin = 450, cmax = 1600, pad = 3.0;   // [m/s], [m/s], [s]
  }

  public static final class Result {
    public double[] asym;      // per trace, NaN = not computed (reference itself, no reference, empty window)
    public double[] s;         // signed position along the line relative to the virtual source [m] (NaN if unknown)
    public double[] dist;      // distance to the virtual source [m]
    public int[] group;        // index of the virtual source (pairing value) of each trace
    public String[] groupNames;
    public double azimuth;     // azimuth of the +s direction [deg from north], NaN if no coordinates
    public boolean hasXY;
  }

  /**
   * @param corr  two-sided correlations [ntr][2nl+1], lag 0 at sample nl
   * @param refOf reference trace of each trace (-1 = none)
   * @param x,y   receiver coordinates (may be null)
   */
  public static Result compute( float[][] corr, int nl, double dt, int[] refOf, double[] x, double[] y, Params p, String[] groupNames, int[] group ) {
    int ntr = corr.length;
    Result r = new Result();
    r.asym = new double[ntr]; r.s = new double[ntr]; r.dist = new double[ntr];
    Arrays.fill( r.asym, Double.NaN ); Arrays.fill( r.s, Double.NaN ); Arrays.fill( r.dist, Double.NaN );
    r.group = group; r.groupNames = groupNames;
    r.hasXY = ( x != null && y != null );
    r.azimuth = Double.NaN;

    // principal direction of the line (all traces with coordinates)
    double ux = 0, uy = 0;
    if( r.hasXY ) {
      double mx = 0, my = 0; int n = 0;
      for( int i = 0; i < ntr; i++ ) { mx += x[i]; my += y[i]; n++; }
      mx /= n; my /= n;
      double sxx = 0, syy = 0, sxy = 0;
      for( int i = 0; i < ntr; i++ ) { double a = x[i]-mx, b = y[i]-my; sxx += a*a; syy += b*b; sxy += a*b; }
      double th = 0.5 * Math.atan2( 2*sxy, sxx - syy );
      ux = Math.cos( th ); uy = Math.sin( th );
      r.azimuth = ( Math.toDegrees( Math.atan2( ux, uy ) ) + 360.0 ) % 360.0;
    }
    for( int i = 0; i < ntr; i++ ) {
      int k = refOf[i];
      if( k < 0 || k == i ) continue;
      int k0 = 1, k1 = nl;
      if( r.hasXY ) {
        double dx = x[i] - x[k], dy = y[i] - y[k];
        r.dist[i] = Math.hypot( dx, dy );
        r.s[i] = dx * ux + dy * uy;
        k0 = Math.max( 1, (int)Math.floor( ( r.dist[i] / p.cmax - p.pad ) / dt ) );
        k1 = Math.min( nl, (int)Math.ceil( ( r.dist[i] / p.cmin + p.pad ) / dt ) );
      }
      else {
        r.s[i] = i - k;      // trace order as a proxy
      }
      if( k1 < k0 ) continue;
      double ep = 0, em = 0;
      float[] c = corr[i];
      for( int j = k0; j <= k1; j++ ) { ep += (double)c[nl + j] * c[nl + j]; em += (double)c[nl - j] * c[nl - j]; }
      if( ep + em > 0 ) r.asym[i] = ( ep - em ) / ( ep + em );
    }
    return r;
  }

  static double median( List<Double> v ) {
    if( v.isEmpty() ) return Double.NaN;
    double[] a = new double[v.size()];
    for( int i = 0; i < a.length; i++ ) a[i] = v.get( i );
    Arrays.sort( a );
    return ( a.length % 2 == 1 ) ? a[a.length/2] : 0.5 * ( a[a.length/2 - 1] + a[a.length/2] );
  }

  /** Text summary (HTML) */
  public static String summary( Result r, Params p ) {
    List<Double> a = new ArrayList<Double>(), d = new ArrayList<Double>();
    int strong = 0;
    for( int i = 0; i < r.asym.length; i++ ) {
      if( Double.isNaN( r.asym[i] ) ) continue;
      a.add( r.asym[i] );
      if( Math.abs( r.asym[i] ) > 0.3 ) strong++;
      if( !Double.isNaN( r.s[i] ) && r.s[i] != 0 ) d.add( r.asym[i] * Math.signum( r.s[i] ) );
    }
    if( a.isEmpty() ) return "Nenhum traço com assimetria calculada.";
    double medA = median( a ), medD = median( d );
    StringBuilder sb = new StringBuilder( "<html>" );
    sb.append( String.format( Locale.US, "<b>%d traços</b> &nbsp; A mediana = <b>%+.2f</b> &nbsp; |A| &gt; 0,3 em <b>%.0f%%</b> &nbsp; D = A·sinal(s) mediana = <b>%+.2f</b><br>",
        a.size(), medA, 100.0 * strong / a.size(), medD ) );
    sb.append( r.hasXY ? String.format( Locale.US, "Janela de chegada: x/%.0f − %.1f s até x/%.0f + %.1f s. ", p.cmax, p.pad, p.cmin, p.pad )
                       : "Sem rec_x/rec_y: janela = lado inteiro, posição = ordem dos traços. " );
    double frac = 100.0 * strong / a.size();
    if( frac < 35 && Math.abs( medD ) < 0.25 ) {
      sb.append( "<font color='#2f6b3a'><b>Campo aproximadamente difuso</b>: somar os dois lados é válido.</font>" );
    }
    else if( !Double.isNaN( medD ) && Math.abs( medD ) >= 0.25 ) {
      if( r.hasXY ) {
        double prop = medD > 0 ? r.azimuth : ( r.azimuth + 180 ) % 360;
        double from = ( prop + 180 ) % 360;
        sb.append( String.format( Locale.US, "<font color='#9b2f2f'><b>Campo direcional</b>: a energia anda ao longo da linha no sentido do azimute %.0f°, ou seja, <b>chega de ~%.0f°</b> (%s). Use só o lado com energia; não some os lados.</font>",
            prop, from, rumo( from ) ) );
      }
      else {
        sb.append( String.format( "<font color='#9b2f2f'><b>Campo direcional</b>: a energia anda no sentido de %s índice de traço. Use só o lado com energia.</font>", medD > 0 ? "maior" : "menor" ) );
      }
    }
    else {
      sb.append( "<font color='#a8661b'><b>Assimetria forte, mas sem sentido único ao longo da linha</b> (fonte próxima ou dentro da linha?). Confira o gráfico.</font>" );
    }
    sb.append( "</html>" );
    return sb.toString();
  }

  static String rumo( double az ) {
    String[] r = { "N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW" };
    return r[(int)Math.round( az / 22.5 ) % 16];
  }

  //====================================================================
  /** Window: A versus signed position along the line, summary, CSV export */
  public static final class Frame extends JFrame {
    private final Result myRes;
    private final JLabel myCursor = new JLabel( " " );
    private static final Color[] COLS = { new Color(15,95,122), new Color(168,102,27), new Color(47,107,58), new Color(155,47,47) };

    public Frame( Result res, Params p, String title, double[] traceNumbers ) {
      super( "Assimetria causal × acausal — " + title );
      myRes = res;
      setDefaultCloseOperation( DISPOSE_ON_CLOSE );
      JLabel info = new JLabel( summary( res, p ) );
      info.setBorder( BorderFactory.createEmptyBorder( 6, 10, 6, 10 ) );
      JButton csv = new JButton( "Exportar CSV..." );
      csv.addActionListener( e -> saveCsv( traceNumbers ) );
      JPanel bar = new JPanel( new FlowLayout( FlowLayout.LEFT, 8, 2 ) );
      JButton help = new JButton( "Ajuda" );
      help.addActionListener( e -> csClockHelp.showInterferometria( this ) );
      bar.add( help ); bar.add( csv ); bar.add( new JLabel( "A > 0: energia indo da fonte virtual para o traço;  A < 0: do traço para a fonte virtual  |" ) ); bar.add( myCursor );
      JPanel top = new JPanel( new BorderLayout() );
      top.add( info, BorderLayout.CENTER ); top.add( bar, BorderLayout.SOUTH );
      Plot plot = new Plot();
      plot.setPreferredSize( new Dimension( 860, 420 ) );
      getContentPane().setLayout( new BorderLayout() );
      getContentPane().add( top, BorderLayout.NORTH );
      getContentPane().add( plot, BorderLayout.CENTER );
      pack();
    }

    private void saveCsv( double[] trc ) {
      JFileChooser fc = new JFileChooser();
      fc.setFileFilter( new FileNameExtensionFilter( "CSV (*.csv)", "csv" ) );
      fc.setSelectedFile( new File( "assimetria.csv" ) );
      if( fc.showSaveDialog( this ) != JFileChooser.APPROVE_OPTION ) return;
      File f = fc.getSelectedFile();
      if( !f.getName().toLowerCase().endsWith( ".csv" ) ) f = new File( f.getPath() + ".csv" );
      try( PrintWriter w = new PrintWriter( f, "UTF-8" ) ) {
        w.println( "traco,fonte_virtual,distancia_m,posicao_s_m,A,D=A*sinal(s)" );
        for( int i = 0; i < myRes.asym.length; i++ ) {
          if( Double.isNaN( myRes.asym[i] ) ) continue;
          String g = ( myRes.groupNames != null && myRes.group[i] >= 0 ) ? myRes.groupNames[myRes.group[i]] : "";
          w.println( String.format( Locale.US, "%.0f,%s,%.1f,%.1f,%.4f,%.4f", trc[i], g, myRes.dist[i], myRes.s[i], myRes.asym[i],
              myRes.asym[i] * Math.signum( myRes.s[i] ) ) );
        }
        myCursor.setText( "Salvo: " + f.getName() );
      }
      catch( IOException ex ) {
        JOptionPane.showMessageDialog( this, "Erro ao salvar CSV:\n" + ex.getMessage() );
      }
    }

    final class Plot extends JPanel {
      static final int ML = 64, MR = 20, MT = 16, MB = 44;
      double smin, smax;
      Plot() {
        setBackground( Color.white );
        smin = Double.POSITIVE_INFINITY; smax = Double.NEGATIVE_INFINITY;
        for( int i = 0; i < myRes.s.length; i++ ) if( !Double.isNaN( myRes.asym[i] ) ) { smin = Math.min( smin, myRes.s[i] ); smax = Math.max( smax, myRes.s[i] ); }
        if( !( smax > smin ) ) { smin = -1; smax = 1; }
        double pad = 0.05 * ( smax - smin ); smin -= pad; smax += pad;
        addMouseMotionListener( new MouseAdapter() {
          @Override public void mouseMoved( MouseEvent e ) {
            int best = -1; double bd = 64;
            for( int i = 0; i < myRes.asym.length; i++ ) {
              if( Double.isNaN( myRes.asym[i] ) ) continue;
              double dx = xp( myRes.s[i] ) - e.getX(), dy = yp( myRes.asym[i] ) - e.getY();
              if( dx*dx + dy*dy < bd ) { bd = dx*dx + dy*dy; best = i; }
            }
            myCursor.setText( best < 0 ? " " : String.format( Locale.US, "traço %d: s = %.0f m, distância %.0f m, A = %+.2f", best + 1, myRes.s[best], myRes.dist[best], myRes.asym[best] ) );
          }
        } );
      }
      double sc() { return myRes.hasXY ? 1000.0 : 1.0; }
      int xp( double s ) { return ML + (int)Math.round( ( s - smin ) / ( smax - smin ) * ( getWidth() - ML - MR ) ); }
      int yp( double a ) { return MT + (int)Math.round( ( 1 - a ) / 2.0 * ( getHeight() - MT - MB ) ); }
      @Override
      protected void paintComponent( Graphics g0 ) {
        super.paintComponent( g0 );
        Graphics2D g = (Graphics2D)g0;
        g.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
        g.setFont( g.getFont().deriveFont( Font.PLAIN, 11f ) );
        FontMetrics fm = g.getFontMetrics();
        int w = getWidth() - ML - MR, h = getHeight() - MT - MB;
        g.setColor( new Color( 250, 250, 250 ) ); g.fillRect( ML, MT, w, h );
        g.setColor( new Color( 225, 235, 240 ) ); g.fillRect( ML, yp( 0.3 ), w, yp( -0.3 ) - yp( 0.3 ) );
        g.setColor( Color.gray );
        for( double a : new double[]{ -1, -0.5, 0, 0.5, 1 } ) {
          int y = yp( a ); g.drawLine( ML - 4, y, ML, y );
          String t = String.format( Locale.US, "%+.1f", a ); g.drawString( t, ML - 8 - fm.stringWidth( t ), y + 4 );
        }
        g.setColor( Color.darkGray ); g.drawLine( ML, yp( 0 ), ML + w, yp( 0 ) );
        if( smin < 0 && smax > 0 ) { g.setColor( Color.gray ); g.drawLine( xp( 0 ), MT, xp( 0 ), MT + h ); g.drawString( "fonte virtual", xp( 0 ) + 4, MT + 12 ); }
        double[] tk = csDispersionFC.ticks( smin / sc(), smax / sc(), Math.max( 2, w / 90 ) );
        for( double t : tk ) {
          int x = xp( t * sc() ); g.setColor( Color.gray ); g.drawLine( x, MT + h, x, MT + h + 4 );
          String s = csDispersionFC.fmt( t ); g.drawString( s, x - fm.stringWidth( s ) / 2, MT + h + 16 );
        }
        String xl = myRes.hasXY ? String.format( Locale.US, "posição ao longo da linha s (km), +s = azimute %.0f°", myRes.azimuth ) : "posição (diferença de índice de traço)";
        g.setColor( Color.darkGray ); g.drawString( xl, ML + w / 2 - fm.stringWidth( xl ) / 2, MT + h + 34 );
        Graphics2D gr = (Graphics2D)g.create(); gr.rotate( -Math.PI / 2 );
        String yl = "A = (E+ − E−)/(E+ + E−)"; gr.drawString( yl, -( MT + h / 2 ) - fm.stringWidth( yl ) / 2, 16 ); gr.dispose();
        g.setStroke( new BasicStroke( 1f ) );
        for( int i = 0; i < myRes.asym.length; i++ ) {
          if( Double.isNaN( myRes.asym[i] ) ) continue;
          Color c = COLS[Math.max( 0, myRes.group[i] ) % COLS.length];
          g.setColor( c ); g.fillOval( xp( myRes.s[i] ) - 4, yp( myRes.asym[i] ) - 4, 8, 8 );
        }
        g.setColor( Color.darkGray ); g.drawRect( ML, MT, w, h );
      }
    }
  }
}
