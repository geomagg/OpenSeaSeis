package cseis.plugins;

import cseis.seaview.csSeisPaneBundle;
import cseis.seaview.plugin.csPluginContext;
import cseis.seaview.plugin.csSeaViewPlugin;
import cseis.seis.csHeader;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicTraceBuffer;
import cseis.seisdisp.csSampleInfo;
import cseis.seisdisp.csSeisView;

import javax.imageio.ImageIO;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingUtilities;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.AWTEvent;
import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.event.AWTEventListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Locale;
import java.util.TimeZone;

/**
 * Espectrograma (STFT) do traço clicado com o mouse.
 * <p>
 * Menu Plugins → "Espectrograma: clique num traço" liga/desliga o modo. Com o modo ligado,
 * um clique com o botão esquerdo em qualquer painel sísmico abre (ou atualiza) a janela
 * do espectrograma daquele traço.
 * <p>
 * Cálculo: remove a média do traço; janelas de Hann com sobreposição; densidade espectral
 * de potência unilateral em dB (10·log10), normalizada pelo máximo do traço (0 dB = máximo).
 */
public class csSpectrogramPlugin implements csSeaViewPlugin {
  private csPluginContext myContext;
  private JCheckBoxMenuItem myToggle;
  private SpectrogramFrame myLastFrame;
  private boolean myNewWindowPerClick = false;
  private final AWTEventListener myClickListener = this::onAwtEvent;

  @Override
  public String getName() { return "Espectrograma"; }

  @Override
  public void install( csPluginContext ctx ) {
    myContext = ctx;
    myToggle = new JCheckBoxMenuItem( "Espectrograma: clique num traço" );
    myToggle.setAccelerator( KeyStroke.getKeyStroke( "ctrl shift S" ) );
    myToggle.addActionListener( e -> setEnabled( myToggle.isSelected() ) );
    ctx.getPluginMenu().add( myToggle );
  }

  private void setEnabled( boolean on ) {
    Toolkit tk = Toolkit.getDefaultToolkit();
    tk.removeAWTEventListener( myClickListener );
    if( on ) {
      tk.addAWTEventListener( myClickListener, AWTEvent.MOUSE_EVENT_MASK );
      myContext.setStatus( "Espectrograma ligado: clique num traço (Ctrl+Shift+S desliga)" );
    }
    else {
      myContext.setStatus( "Espectrograma desligado" );
    }
    if( myToggle.isSelected() != on ) myToggle.setSelected( on );
  }

  private void onAwtEvent( AWTEvent event ) {
    if( event.getID() != MouseEvent.MOUSE_CLICKED ) return;
    MouseEvent me = (MouseEvent)event;
    if( !SwingUtilities.isLeftMouseButton( me ) || !( me.getSource() instanceof csSeisView ) ) return;
    csSeisView view = (csSeisView)me.getSource();
    csISeismicTraceBuffer buf = view.getTraceBuffer();
    if( buf == null || buf.numTraces() == 0 ) return;
    csSampleInfo info = view.getSampleInfo( me.getX(), me.getY() );
    int itrc = info.trace;
    if( itrc < 0 || itrc >= buf.numTraces() ) return;

    csSeisPaneBundle bundle = (csSeisPaneBundle)SwingUtilities.getAncestorOfClass( csSeisPaneBundle.class, view );
    csHeaderDef[] defs = bundle != null ? bundle.getTraceHeaderDefs() : new csHeaderDef[0];
    String fileName = bundle != null ? new File( bundle.getFilenamePath() ).getName() : "";
    float[] s = buf.samples( itrc );
    TraceData td = new TraceData( s.clone(), view.getSampleInt() / 1000.0, itrc, buf.originalTraceNumber( itrc ),
                                  describeTrace( buf.headerValues( itrc ), defs ), startTimeUtc( buf.headerValues( itrc ), defs ), fileName );

    if( myLastFrame == null || !myLastFrame.isDisplayable() || myNewWindowPerClick ) {
      SpectrogramFrame f = new SpectrogramFrame( myLastFrame );
      if( myLastFrame != null && myLastFrame.isDisplayable() ) {
        f.setLocation( myLastFrame.getX() + 30, myLastFrame.getY() + 30 );
      }
      else {
        f.setLocationRelativeTo( myContext.getSeaView() );
      }
      myLastFrame = f;
    }
    myLastFrame.setTrace( td );
    myLastFrame.setVisible( true );
    myLastFrame.toFront();
  }

  //------------------------------------------------------------------
  private static final String[] PREFERRED_HEADERS = { "trcno", "fileno", "rec_line", "rcv", "node", "chan", "offset", "cdp" };

  private static String describeTrace( csHeader[] h, csHeaderDef[] defs ) {
    StringBuilder sb = new StringBuilder();
    if( h == null ) return "";
    for( String name : PREFERRED_HEADERS ) {
      for( int k = 0; k < defs.length && k < h.length; k++ ) {
        if( defs[k].name.equalsIgnoreCase( name ) && h[k] != null ) {
          if( sb.length() > 0 ) sb.append( "  |  " );
          sb.append( name ).append( ' ' ).append( h[k].value() );
        }
      }
    }
    return sb.toString();
  }

  /** Início do traço em UTC a partir de time_samp1 (+ time_samp1_us), ou null. */
  private static String startTimeUtc( csHeader[] h, csHeaderDef[] defs ) {
    if( h == null ) return null;
    double t = -1, us = 0;
    for( int k = 0; k < defs.length && k < h.length; k++ ) {
      if( h[k] == null ) continue;
      if( defs[k].name.equalsIgnoreCase( "time_samp1" ) ) t = h[k].doubleValue();
      else if( defs[k].name.equalsIgnoreCase( "time_samp1_us" ) ) us = h[k].doubleValue();
    }
    if( t <= 0 ) return null;
    SimpleDateFormat f = new SimpleDateFormat( "yyyy-MM-dd HH:mm:ss.SSS 'UTC'", Locale.US );
    f.setTimeZone( TimeZone.getTimeZone( "UTC" ) );
    return f.format( new Date( (long)( t * 1000.0 + us / 1000.0 ) ) );
  }

  //==================================================================
  static final class TraceData {
    final float[] samples;
    final double dt;            // [s]
    final int traceIndex;
    final int traceNumber;
    final String headers;
    final String startUtc;
    final String fileName;
    TraceData( float[] samples, double dt, int traceIndex, int traceNumber, String headers, String startUtc, String fileName ) {
      this.samples = samples; this.dt = dt; this.traceIndex = traceIndex; this.traceNumber = traceNumber;
      this.headers = headers; this.startUtc = startUtc; this.fileName = fileName;
    }
  }

  /** Resultado da STFT: psd[frame][bin] em dB (não normalizado). */
  static final class Stft {
    double[][] psd;
    double[] frameTime;   // centro de cada janela [s]
    double df;            // [Hz]
    double maxDb;
    int nfreq;
  }

  static Stft compute( float[] x, double dt, int nwin, int overlapPct ) {
    int n = x.length;
    if( nwin > n ) nwin = Integer.highestOneBit( Math.max( n, 2 ) );
    int hop = Math.max( 1, (int)Math.round( nwin * ( 1.0 - overlapPct / 100.0 ) ) );
    int nframes = Math.max( 1, 1 + ( n - nwin ) / hop );
    int nfreq = nwin / 2 + 1;

    double mean = 0;
    for( float v : x ) mean += v;
    mean /= Math.max( 1, n );

    double[] w = new double[nwin];
    double wss = 0;
    for( int i = 0; i < nwin; i++ ) {
      w[i] = 0.5 - 0.5 * Math.cos( 2.0 * Math.PI * i / nwin );
      wss += w[i] * w[i];
    }
    double fs = 1.0 / dt;
    double scale = 1.0 / ( fs * wss );

    Stft r = new Stft();
    r.psd = new double[nframes][nfreq];
    r.frameTime = new double[nframes];
    r.df = fs / nwin;
    r.nfreq = nfreq;
    r.maxDb = -Double.MAX_VALUE;
    double[] re = new double[nwin];
    double[] im = new double[nwin];
    for( int f = 0; f < nframes; f++ ) {
      int i0 = f * hop;
      for( int i = 0; i < nwin; i++ ) {
        int k = i0 + i;
        re[i] = k < n ? ( x[k] - mean ) * w[i] : 0.0;
        im[i] = 0.0;
      }
      fft( re, im );
      for( int b = 0; b < nfreq; b++ ) {
        double p = ( re[b] * re[b] + im[b] * im[b] ) * scale;
        if( b > 0 && b < nwin / 2 ) p *= 2.0;          // unilateral
        double db = 10.0 * Math.log10( Math.max( p, 1.0e-30 ) );
        r.psd[f][b] = db;
        if( db > r.maxDb ) r.maxDb = db;
      }
      r.frameTime[f] = ( i0 + nwin / 2.0 ) * dt;
    }
    return r;
  }

  /** FFT radix-2 in-place (n potência de 2). */
  static void fft( double[] re, double[] im ) {
    int n = re.length;
    for( int i = 1, j = 0; i < n; i++ ) {
      int bit = n >> 1;
      for( ; ( j & bit ) != 0; bit >>= 1 ) j ^= bit;
      j ^= bit;
      if( i < j ) {
        double t = re[i]; re[i] = re[j]; re[j] = t;
        t = im[i]; im[i] = im[j]; im[j] = t;
      }
    }
    for( int len = 2; len <= n; len <<= 1 ) {
      double ang = -2.0 * Math.PI / len;
      double wr = Math.cos( ang ), wi = Math.sin( ang );
      for( int i = 0; i < n; i += len ) {
        double cr = 1.0, ci = 0.0;
        for( int k = 0; k < len / 2; k++ ) {
          int a = i + k, b = i + k + len / 2;
          double tr = re[b] * cr - im[b] * ci;
          double ti = re[b] * ci + im[b] * cr;
          re[b] = re[a] - tr; im[b] = im[a] - ti;
          re[a] += tr;        im[a] += ti;
          double ncr = cr * wr - ci * wi;
          ci = cr * wi + ci * wr;
          cr = ncr;
        }
      }
    }
  }

  //------------------------------------------------------------------
  // Paleta tipo "viridis" (interpolação linear entre pontos de controle)
  private static final float[][] VIRIDIS = {
    { 0.267f, 0.005f, 0.329f }, { 0.283f, 0.141f, 0.458f }, { 0.254f, 0.265f, 0.530f },
    { 0.207f, 0.372f, 0.553f }, { 0.164f, 0.471f, 0.558f }, { 0.128f, 0.567f, 0.551f },
    { 0.135f, 0.659f, 0.518f }, { 0.267f, 0.749f, 0.441f }, { 0.478f, 0.821f, 0.318f },
    { 0.741f, 0.873f, 0.150f }, { 0.993f, 0.906f, 0.144f } };

  static int colorRgb( double v ) {        // v em [0,1]
    if( !( v > 0 ) ) v = 0;
    if( v > 1 ) v = 1;
    double p = v * ( VIRIDIS.length - 1 );
    int i = Math.min( (int)p, VIRIDIS.length - 2 );
    double t = p - i;
    int r = (int)( 255 * ( VIRIDIS[i][0] + t * ( VIRIDIS[i+1][0] - VIRIDIS[i][0] ) ) );
    int g = (int)( 255 * ( VIRIDIS[i][1] + t * ( VIRIDIS[i+1][1] - VIRIDIS[i][1] ) ) );
    int b = (int)( 255 * ( VIRIDIS[i][2] + t * ( VIRIDIS[i+1][2] - VIRIDIS[i][2] ) ) );
    return ( r << 16 ) | ( g << 8 ) | b;
  }

  //==================================================================
  final class SpectrogramFrame extends JFrame {
    private static final int ML = 70, MR = 90, MT = 12, MB = 44;   // margens do gráfico
    private TraceData myTrace;
    private Stft myStft;
    private final PlotPanel myPlot = new PlotPanel();
    private final JLabel myInfo = new JLabel( " " );
    private final JLabel myCursor = new JLabel( " " );
    private final JComboBox<Integer> myWin;
    private final JSpinner myOverlap = new JSpinner( new SpinnerNumberModel( 75, 0, 95, 5 ) );
    private final JTextField myFmin = new JTextField( "0", 5 );
    private final JTextField myFmax = new JTextField( "", 5 );
    private final JComboBox<String> myFScale = new JComboBox<>( new String[]{ "linear", "log" } );
    private final JComboBox<Integer> myRange = new JComboBox<>( new Integer[]{ 40, 60, 80, 100, 120 } );

    SpectrogramFrame( SpectrogramFrame copySettingsFrom ) {
      super( "Espectrograma" );
      setDefaultCloseOperation( DISPOSE_ON_CLOSE );
      List<Integer> sizes = new ArrayList<>();
      for( int s = 32; s <= 16384; s *= 2 ) sizes.add( s );
      myWin = new JComboBox<>( sizes.toArray( new Integer[0] ) );
      myRange.setSelectedItem( 60 );
      if( copySettingsFrom != null ) {
        myWin.setSelectedItem( copySettingsFrom.myWin.getSelectedItem() );
        myOverlap.setValue( copySettingsFrom.myOverlap.getValue() );
        myFmin.setText( copySettingsFrom.myFmin.getText() );
        myFmax.setText( copySettingsFrom.myFmax.getText() );
        myFScale.setSelectedItem( copySettingsFrom.myFScale.getSelectedItem() );
        myRange.setSelectedItem( copySettingsFrom.myRange.getSelectedItem() );
      }

      JPanel controls = new JPanel( new FlowLayout( FlowLayout.LEFT, 6, 2 ) );
      controls.add( new JLabel( "Janela (amostras):" ) ); controls.add( myWin );
      controls.add( new JLabel( "Sobreposição %:" ) );   controls.add( myOverlap );
      controls.add( new JLabel( "f mín:" ) ); controls.add( myFmin );
      controls.add( new JLabel( "f máx [Hz]:" ) ); controls.add( myFmax );
      controls.add( new JLabel( "Escala f:" ) ); controls.add( myFScale );
      controls.add( new JLabel( "Faixa [dB]:" ) ); controls.add( myRange );

      JPanel controls2 = new JPanel( new FlowLayout( FlowLayout.LEFT, 6, 2 ) );
      JCheckBox newWin = new JCheckBox( "Nova janela a cada clique (comparar traços)", myNewWindowPerClick );
      newWin.addActionListener( e -> myNewWindowPerClick = newWin.isSelected() );
      JButton save = new JButton( "Salvar PNG..." );
      save.addActionListener( e -> savePng() );
      controls2.add( newWin );
      controls2.add( save );
      controls2.add( myCursor );

      JPanel top = new JPanel( new GridLayout( 3, 1 ) );
      myInfo.setBorder( BorderFactory.createEmptyBorder( 2, 8, 2, 8 ) );
      top.add( myInfo ); top.add( controls ); top.add( controls2 );

      getContentPane().setLayout( new BorderLayout() );
      getContentPane().add( top, BorderLayout.NORTH );
      getContentPane().add( myPlot, BorderLayout.CENTER );
      myPlot.setPreferredSize( new Dimension( 900, 560 ) );

      myWin.addActionListener( e -> recompute() );
      myOverlap.addChangeListener( e -> recompute() );
      myFmin.addActionListener( e -> myPlot.repaint() );
      myFmax.addActionListener( e -> myPlot.repaint() );
      myFScale.addActionListener( e -> myPlot.repaint() );
      myRange.addActionListener( e -> myPlot.repaint() );
      pack();
    }

    void setTrace( TraceData td ) {
      boolean first = myTrace == null;
      myTrace = td;
      if( first && myFmax.getText().trim().isEmpty() ) {
        myFmax.setText( fmt( 0.5 / td.dt ) );
        int ns = td.samples.length;
        int w = Integer.highestOneBit( Math.max( 32, ns / 16 ) );
        myWin.setSelectedItem( Math.min( 4096, w ) );    // dispara recompute()
      }
      setTitle( "Espectrograma — traço " + td.traceNumber + ( td.fileName.isEmpty() ? "" : "  (" + td.fileName + ")" ) );
      recompute();
    }

    private void recompute() {
      if( myTrace == null ) return;
      int nwin = (Integer)myWin.getSelectedItem();
      int ov = (Integer)myOverlap.getValue();
      myStft = compute( myTrace.samples, myTrace.dt, nwin, ov );
      int nwinUsed = 2 * ( myStft.nfreq - 1 );
      myInfo.setText( String.format( Locale.US, "<html><b>Traço %d</b> &nbsp;&nbsp; %s%s<br>%d amostras, dt = %s ms, duração %s s &nbsp;&nbsp;|&nbsp;&nbsp; janela %d amostras = %s s, Δf = %s Hz, %d janelas</html>",
          myTrace.traceNumber, myTrace.headers.replace( "  |  ", " &nbsp;|&nbsp; " ), myTrace.startUtc == null ? "" : " &nbsp;&nbsp; início: " + myTrace.startUtc,
          myTrace.samples.length, fmt( myTrace.dt * 1000 ), fmt( myTrace.samples.length * myTrace.dt ),
          nwinUsed, fmt( nwinUsed * myTrace.dt ), fmt( myStft.df ), myStft.psd.length ) );
      myPlot.invalidateImage();
      myPlot.repaint();
    }

    double fNyq() { return 0.5 / myTrace.dt; }

    double fMin() {
      double f = parse( myFmin.getText(), 0 );
      if( isLog() && f <= 0 ) f = myStft.df;       // log: começa no primeiro bin
      return Math.max( 0, Math.min( f, fNyq() ) );
    }
    double fMax() {
      double f = parse( myFmax.getText(), fNyq() );
      if( f <= fMin() ) f = fNyq();
      return Math.min( f, fNyq() );
    }
    boolean isLog() { return "log".equals( myFScale.getSelectedItem() ); }
    int rangeDb() { return (Integer)myRange.getSelectedItem(); }

    /** Frequência na linha y (0 = topo) de um gráfico de altura h. */
    double freqAt( double yFrac ) {                  // yFrac: 0 = base (fmin), 1 = topo (fmax)
      double a = fMin(), b = fMax();
      if( isLog() ) return Math.exp( Math.log( a ) + yFrac * ( Math.log( b ) - Math.log( a ) ) );
      return a + yFrac * ( b - a );
    }
    double yFracOf( double f ) {
      double a = fMin(), b = fMax();
      if( isLog() ) return ( Math.log( f ) - Math.log( a ) ) / ( Math.log( b ) - Math.log( a ) );
      return ( f - a ) / ( b - a );
    }

    private void savePng() {
      if( myTrace == null ) return;
      JFileChooser fc = new JFileChooser();
      fc.setFileFilter( new FileNameExtensionFilter( "PNG (*.png)", "png" ) );
      fc.setSelectedFile( new File( "espectrograma_traco" + myTrace.traceNumber + ".png" ) );
      if( fc.showSaveDialog( this ) != JFileChooser.APPROVE_OPTION ) return;
      File f = fc.getSelectedFile();
      if( !f.getName().toLowerCase().endsWith( ".png" ) ) f = new File( f.getPath() + ".png" );
      BufferedImage img = new BufferedImage( myPlot.getWidth(), myPlot.getHeight(), BufferedImage.TYPE_INT_RGB );
      Graphics2D g = img.createGraphics();
      myPlot.paint( g );
      g.dispose();
      try {
        ImageIO.write( img, "png", f );
        myCursor.setText( "Salvo: " + f.getName() );
      }
      catch( IOException ex ) {
        myContext.error( "Erro ao salvar PNG:\n" + ex.getMessage() );
      }
    }

    //----------------------------------------------------------------
    final class PlotPanel extends JPanel {
      private BufferedImage myImage;
      private String myImageKey = "";

      PlotPanel() {
        setBackground( Color.white );
        MouseAdapter ma = new MouseAdapter() {
          @Override public void mouseMoved( MouseEvent e ) { updateCursor( e.getX(), e.getY() ); }
          @Override public void mouseExited( MouseEvent e ) { myCursor.setText( " " ); }
        };
        addMouseMotionListener( ma );
        addMouseListener( ma );
      }

      void invalidateImage() { myImageKey = ""; }

      // Áreas: forma de onda em cima (~22%), espectrograma embaixo
      int waveTop()    { return MT; }
      int waveHeight() { return Math.max( 60, (int)( ( getHeight() - MT - MB ) * 0.22 ) ); }
      int specTop()    { return waveTop() + waveHeight() + 14; }
      int specHeight() { return Math.max( 20, getHeight() - MB - specTop() ); }
      int plotWidth()  { return Math.max( 20, getWidth() - ML - MR ); }
      double duration() { return myTrace.samples.length * myTrace.dt; }

      private void updateCursor( int x, int y ) {
        if( myStft == null ) return;
        int px = x - ML, py = y - specTop();
        if( px < 0 || px >= plotWidth() ) { myCursor.setText( " " ); return; }
        double t = px * duration() / plotWidth();
        if( py >= 0 && py < specHeight() ) {
          double f = freqAt( 1.0 - ( py + 0.5 ) / specHeight() );
          double db = valueAt( t, f ) - myStft.maxDb;
          myCursor.setText( String.format( Locale.US, "t = %s s   f = %s Hz   %.1f dB", fmt( t ), fmt( f ), db ) );
        }
        else {
          int is = Math.min( myTrace.samples.length - 1, (int)( t / myTrace.dt ) );
          myCursor.setText( String.format( Locale.US, "t = %s s   amplitude = %g", fmt( t ), myTrace.samples[is] ) );
        }
      }

      /** Valor (dB) da STFT no tempo t e frequência f: janela mais próxima, interpolação linear em f. */
      double valueAt( double t, double f ) {
        double[] ft = myStft.frameTime;
        int fr;
        if( ft.length == 1 ) fr = 0;
        else {
          double step = ft[1] - ft[0];
          fr = (int)Math.round( ( t - ft[0] ) / step );
          fr = Math.max( 0, Math.min( ft.length - 1, fr ) );
        }
        double bin = f / myStft.df;
        int b0 = (int)Math.floor( bin );
        if( b0 >= myStft.nfreq - 1 ) return myStft.psd[fr][myStft.nfreq - 1];
        if( b0 < 0 ) return myStft.psd[fr][0];
        double w = bin - b0;
        return ( 1 - w ) * myStft.psd[fr][b0] + w * myStft.psd[fr][b0 + 1];
      }

      private void buildImage( int w, int h ) {
        String key = w + "x" + h + "|" + fMin() + "|" + fMax() + "|" + isLog() + "|" + rangeDb();
        if( myImage != null && key.equals( myImageKey ) ) return;
        myImage = new BufferedImage( w, h, BufferedImage.TYPE_INT_RGB );
        double dur = duration();
        double top = myStft.maxDb, range = rangeDb();
        double[] freqs = new double[h];
        for( int y = 0; y < h; y++ ) freqs[y] = freqAt( 1.0 - ( y + 0.5 ) / h );
        for( int x = 0; x < w; x++ ) {
          double t = ( x + 0.5 ) * dur / w;
          for( int y = 0; y < h; y++ ) {
            double v = ( valueAt( t, freqs[y] ) - ( top - range ) ) / range;
            myImage.setRGB( x, y, colorRgb( v ) );
          }
        }
        myImageKey = key;
      }

      @Override
      protected void paintComponent( Graphics g0 ) {
        super.paintComponent( g0 );
        Graphics2D g = (Graphics2D)g0;
        g.setColor( Color.white );
        g.fillRect( 0, 0, getWidth(), getHeight() );
        if( myTrace == null || myStft == null ) {
          g.setColor( Color.gray );
          g.drawString( "Clique num traço do SeaView.", ML, MT + 20 );
          return;
        }
        g.setRenderingHint( RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON );
        g.setFont( g.getFont().deriveFont( Font.PLAIN, 11f ) );
        FontMetrics fm = g.getFontMetrics();
        int pw = plotWidth();
        double dur = duration();

        //--- Forma de onda
        int wt = waveTop(), wh = waveHeight();
        float amax = 0;
        for( float v : myTrace.samples ) amax = Math.max( amax, Math.abs( v ) );
        if( amax == 0 ) amax = 1;
        g.setColor( new Color( 245, 245, 245 ) );
        g.fillRect( ML, wt, pw, wh );
        g.setColor( new Color( 40, 40, 40 ) );
        int ns = myTrace.samples.length;
        int yMid = wt + wh / 2;
        if( ns > 2 * pw ) {                         // min/max por coluna de pixel
          for( int x = 0; x < pw; x++ ) {
            int i0 = (int)( (long)x * ns / pw ), i1 = Math.max( i0 + 1, (int)( (long)( x + 1 ) * ns / pw ) );
            float mn = Float.MAX_VALUE, mx = -Float.MAX_VALUE;
            for( int i = i0; i < i1 && i < ns; i++ ) { mn = Math.min( mn, myTrace.samples[i] ); mx = Math.max( mx, myTrace.samples[i] ); }
            g.drawLine( ML + x, yMid - (int)( mx / amax * wh / 2 ), ML + x, yMid - (int)( mn / amax * wh / 2 ) );
          }
        }
        else {
          int px = -1, py = 0;
          for( int i = 0; i < ns; i++ ) {
            int x = ML + (int)( ( i + 0.5 ) * pw / ns );
            int y = yMid - (int)( myTrace.samples[i] / amax * wh / 2 );
            if( px >= 0 ) g.drawLine( px, py, x, y );
            px = x; py = y;
          }
        }
        g.setColor( Color.gray );
        g.drawRect( ML, wt, pw, wh );
        drawRight( g, fm, String.format( Locale.US, "%.3g", amax ), ML - 4, wt + fm.getAscent() );
        drawRight( g, fm, String.format( Locale.US, "%.3g", -amax ), ML - 4, wt + wh );

        //--- Espectrograma
        int st = specTop(), sh = specHeight();
        buildImage( pw, sh );
        g.drawImage( myImage, ML, st, null );
        g.setColor( Color.darkGray );
        g.drawRect( ML, st, pw, sh );

        // Eixo de tempo
        double[] tt = niceTicks( 0, dur, Math.max( 2, pw / 90 ) );
        for( double t : tt ) {
          int x = ML + (int)Math.round( t / dur * pw );
          g.setColor( Color.darkGray );
          g.drawLine( x, st + sh, x, st + sh + 4 );
          g.drawLine( x, wt + wh, x, wt + wh + 3 );
          String s = fmt( t );
          g.drawString( s, x - fm.stringWidth( s ) / 2, st + sh + 4 + fm.getAscent() );
        }
        String tl = "Tempo [s]";
        g.drawString( tl, ML + pw / 2 - fm.stringWidth( tl ) / 2, st + sh + 6 + 2 * fm.getHeight() );

        // Eixo de frequência
        double[] ft = isLog() ? logTicks( fMin(), fMax() ) : niceTicks( fMin(), fMax(), Math.max( 2, sh / 45 ) );
        for( double f : ft ) {
          int y = st + sh - (int)Math.round( yFracOf( f ) * sh );
          if( y < st - 1 || y > st + sh + 1 ) continue;
          g.setColor( Color.darkGray );
          g.drawLine( ML - 4, y, ML, y );
          drawRight( g, fm, fmt( f ), ML - 6, y + fm.getAscent() / 2 - 1 );
        }
        Graphics2D gr = (Graphics2D)g.create();
        gr.rotate( -Math.PI / 2 );
        String fl = "Frequência [Hz]";
        gr.drawString( fl, -( st + sh / 2 ) - fm.stringWidth( fl ) / 2, 14 );
        gr.dispose();

        // Barra de cores
        int cbx = ML + pw + 18, cbw = 16;
        for( int y = 0; y < sh; y++ ) {
          g.setColor( new Color( colorRgb( 1.0 - (double)y / ( sh - 1 ) ) ) );
          g.drawLine( cbx, st + y, cbx + cbw, st + y );
        }
        g.setColor( Color.darkGray );
        g.drawRect( cbx, st, cbw, sh );
        int range = rangeDb();
        int stepDb = range <= 60 ? 10 : 20;
        for( int d = 0; d <= range; d += stepDb ) {
          int y = st + (int)Math.round( (double)d / range * sh );
          g.drawLine( cbx + cbw, y, cbx + cbw + 3, y );
          g.drawString( ( d == 0 ? "0" : "-" + d ), cbx + cbw + 5, y + fm.getAscent() / 2 - 1 );
        }
        g.drawString( "dB", cbx, st - 2 );
        g.setStroke( new BasicStroke( 1f ) );
      }

      private void drawRight( Graphics g, FontMetrics fm, String s, int xRight, int y ) {
        g.drawString( s, xRight - fm.stringWidth( s ), y );
      }
    }
  }

  //------------------------------------------------------------------
  static double[] niceTicks( double a, double b, int maxTicks ) {
    double span = b - a;
    if( !( span > 0 ) ) return new double[]{ a };
    double raw = span / maxTicks;
    double mag = Math.pow( 10, Math.floor( Math.log10( raw ) ) );
    double step = mag;
    for( double m : new double[]{ 1, 2, 2.5, 5, 10 } ) {
      if( m * mag >= raw ) { step = m * mag; break; }
    }
    List<Double> list = new ArrayList<>();
    for( double v = Math.ceil( a / step - 1e-9 ) * step; v <= b + 1e-9 * step; v += step ) list.add( Math.abs( v ) < 1e-12 ? 0 : v );
    double[] r = new double[list.size()];
    for( int i = 0; i < r.length; i++ ) r[i] = list.get( i );
    return r;
  }

  static double[] logTicks( double a, double b ) {
    List<Double> list = new ArrayList<>();
    int e0 = (int)Math.floor( Math.log10( a ) ), e1 = (int)Math.ceil( Math.log10( b ) );
    boolean fewDecades = ( e1 - e0 ) <= 2;
    for( int e = e0; e <= e1; e++ ) {
      for( int m : fewDecades ? new int[]{ 1, 2, 5 } : new int[]{ 1 } ) {
        double v = m * Math.pow( 10, e );
        if( v >= a * 0.999 && v <= b * 1.001 ) list.add( v );
      }
    }
    double[] r = new double[list.size()];
    for( int i = 0; i < r.length; i++ ) r[i] = list.get( i );
    return r;
  }

  static double parse( String s, double def ) {
    try { return Double.parseDouble( s.trim().replace( ',', '.' ) ); }
    catch( Exception e ) { return def; }
  }

  static String fmt( double v ) {
    if( v == Math.rint( v ) && Math.abs( v ) < 1e9 ) return String.valueOf( (long)v );
    String s = String.format( Locale.US, "%.4g", v );
    if( s.contains( "." ) && !s.contains( "e" ) ) s = s.replaceAll( "0+$", "" ).replaceAll( "\\.$", "" );
    return s;
  }

  /** Para testes sem interface: calcula a STFT de um vetor. */
  static Stft computeForTest( float[] x, double dt, int nwin, int ov ) { return compute( x, dt, nwin, ov ); }
}
