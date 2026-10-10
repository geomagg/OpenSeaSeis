/* SeaView plugin Interferometria: virtual shot gather (VSG) by cross-correlation with a reference trace
   M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import cseis.jni.csVirtualSeismicReader;
import cseis.seis.csHeader;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicTraceBuffer;
import cseis.seis.csTraceBuffer;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.util.Arrays;
import java.util.Locale;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 * Seismic interferometry of the displayed traces: every trace of the active pane is cross-correlated
 * with a reference trace (virtual source); the correlations form a virtual shot gather (VSG) shown in a new pane.
 * <p>
 * Per trace and per stacking window: demean, taper, band-pass (cosine-tapered mask fmin-fmax),
 * temporal normalization (none, one-bit, running absolute mean), spectral whitening (optional),
 * then C(tau) = IFFT( conj(R) X ), stacked over the windows. tau &gt; 0: the trace lags the reference.
 * Output: lags -maxlag..+maxlag (lag 0 at sample maxlag), or the sum of both sides C(tau)+C(-tau),
 * tau = 0..maxlag. Headers are copied; optionally 'offset' = distance to the reference (rec_x, rec_y).
 */
public class csPluginInterferometria implements csISeaViewPlugin {
  static final String NONE = "(nenhum)";
  static final String[] NORMS = { "nenhuma", "one-bit", "média absoluta móvel" };
  private static int ourCounter = 0;
  private final Params myParams = new Params();
  private String myRefHdr = null, myRefVal = null, myPairHdr = null;
  private boolean myDisp = false;
  private boolean myAsym = false;
  private final csAssimetria.Params myAsymParams = new csAssimetria.Params();

  /** Computation parameters (all times in seconds, frequencies in Hz) */
  public static final class Params {
    public double maxLag = 2.0;
    public double fmin = 0.0, fmax = 0.0;       // fmax <= 0: Nyquist
    public int    tempNorm = 0;                 // 0 none, 1 one-bit, 2 running absolute mean
    public double ramWin = 0.5;                 // running absolute mean window [s]
    public boolean whiten = false;
    public double winLen = 0.0;                 // stacking window [s], 0 = whole trace
    public double overlap = 0.0;                // overlap between stacking windows [%]
    public boolean symmetric = false;           // sum of causal and acausal sides
    public boolean normalize = true;            // each output trace / max |C|
  }

  @Override
  public String getName() {
    return "Interferometria (virtual shot gather)...";
  }
  @Override
  public String getDescription() {
    return "Virtual shot gather (VSG): correlação cruzada de cada traço exibido com o traço de referência (fonte virtual), com one-bit, branqueamento e soma dos dois lados";
  }

  //====================================================================
  @Override
  public void run( csIPluginContext ctx ) {
    csISeismicTraceBuffer buf = ctx.getTraceBuffer();
    csHeaderDef[] defs = ctx.getHeaderDef();
    if( buf == null || buf.numTraces() == 0 || defs == null ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "Nenhum dado no painel ativo", "Interferometria", JOptionPane.WARNING_MESSAGE );
      return;
    }
    double dt = ctx.getSampleInt() / 1000.0;
    int ns = buf.numSamples();
    double fnyq = 0.5 / dt;

    String[] names = new String[defs.length];
    for( int i = 0; i < defs.length; i++ ) names[i] = defs[i].name;
    String[] names2 = new String[defs.length+1];
    names2[0] = NONE;
    System.arraycopy( names, 0, names2, 1, names.length );

    //--- Defaults
    if( myRefHdr == null || ctx.getHeaderIndex( myRefHdr ) < 0 ) {
      myRefHdr = ( ctx.getHeaderIndex("rcv") >= 0 ) ? "rcv" : ( ctx.getHeaderIndex("trcno") >= 0 ? "trcno" : names[0] );
      myRefVal = null;
    }
    if( myRefVal == null ) myRefVal = buf.headerValues(0)[ ctx.getHeaderIndex( myRefHdr ) ].toString();
    if( myPairHdr == null ) myPairHdr = ( ctx.getHeaderIndex("chan") >= 0 ) ? "chan" : NONE;
    if( myParams.fmax <= 0 || myParams.fmax > fnyq ) myParams.fmax = fnyq;

    JComboBox<String> comboRef  = new JComboBox<String>( names );
    JTextField textRefVal       = new JTextField( myRefVal, 8 );
    JComboBox<String> comboPair = new JComboBox<String>( names2 );
    JComboBox<String> comboNorm = new JComboBox<String>( NORMS );
    JTextField textRam    = new JTextField( fmt( myParams.ramWin ), 6 );
    JCheckBox  boxWhiten  = new JCheckBox( "Branqueamento espectral", myParams.whiten );
    JTextField textFmin   = new JTextField( fmt( myParams.fmin ), 6 );
    JTextField textFmax   = new JTextField( fmt( myParams.fmax ), 6 );
    JTextField textLag    = new JTextField( fmt( Math.min( myParams.maxLag, (ns-1)*dt ) ), 6 );
    JTextField textWin    = new JTextField( fmt( myParams.winLen ), 6 );
    JTextField textOverlap = new JTextField( fmt( myParams.overlap ), 6 );
    JCheckBox  boxSym     = new JCheckBox( "Somar os dois lados: C(t) + C(-t)", myParams.symmetric );
    JCheckBox  boxNormOut = new JCheckBox( "Normalizar cada correlação (máx = 1)", myParams.normalize );
    // offset is created in the VSG pane if the data do not have it (e.g. data read by INPUT_HDF5)
    boolean hasXY = ctx.getHeaderIndex("rec_x") >= 0 && ctx.getHeaderIndex("rec_y") >= 0;
    JCheckBox  boxOffset  = new JCheckBox( "offset = distância à fonte virtual (rec_x, rec_y)", hasXY );
    boxOffset.setEnabled( hasXY );
    JCheckBox  boxDisp    = new JCheckBox( "Imagem de dispersão f-c do VSG (velocidade de fase x frequência)", hasXY && myDisp );
    boxDisp.setEnabled( hasXY );
    boxDisp.setToolTipText( "Phase-shift sobre as correlações do VSG, usando a distância à fonte virtual (rec_x, rec_y). Também disponível no menu: Dispersão f-c (painel ativo)" );
    JCheckBox  boxAsym    = new JCheckBox( "Assimetria causal × acausal  A = (E+ − E−)/(E+ + E−)", myAsym );
    boxAsym.setToolTipText( "Compara a energia dos lags positivos e negativos na janela de chegada (x/c máx − folga ... x/c mín + folga). Campo difuso: A ~ 0. Grava o cabeçalho 'asym' e abre um gráfico de A ao longo da linha" );
    JTextField textAsCmin = new JTextField( fmt( myAsymParams.cmin ), 5 ), textAsCmax = new JTextField( fmt( myAsymParams.cmax ), 5 );
    JTextField textAsPad  = new JTextField( fmt( myAsymParams.pad ), 4 );
    textAsCmin.setToolTipText( "Velocidade mínima da janela de chegada [m/s] (ex.: Scholte 450; modos da água 1500)" );
    textAsCmax.setToolTipText( "Velocidade máxima da janela de chegada [m/s] (ex.: Scholte 1600; modos da água 3000)" );
    textAsPad.setToolTipText( "Folga [s] antes e depois da janela (~ um período)" );
    comboRef.setSelectedItem( myRefHdr );
    comboPair.setSelectedItem( ctx.getHeaderIndex( myPairHdr ) >= 0 ? myPairHdr : NONE );
    comboNorm.setSelectedIndex( myParams.tempNorm );
    comboRef.setMaximumRowCount( 20 );
    comboPair.setMaximumRowCount( 20 );
    comboRef.addActionListener( e -> {
      int ih = ctx.getHeaderIndex( (String)comboRef.getSelectedItem() );
      if( ih >= 0 ) textRefVal.setText( buf.headerValues(0)[ih].toString() );
    });
    textRam.setToolTipText( "Janela da média absoluta móvel [s] (ex.: metade do maior período de interesse)" );
    textOverlap.setToolTipText( "Sobreposição entre janelas consecutivas [%] (0 a 90). Ex.: 5 min, janelas de 60 s, 50% -> 9 janelas. Com a janela quase retangular usada aqui, a sobreposição quase não melhora a razão sinal/ruído e aumenta o tempo" );
    textWin.setToolTipText( "Os traços são cortados em janelas deste comprimento e as correlações empilhadas. 0 = traço inteiro" );
    comboPair.setToolTipText( "Cada traço é correlacionado com a referência que tem o mesmo valor deste cabeçalho (ex.: chan = componente)" );

    JPanel p = new JPanel( new GridBagLayout() );
    int row = 0;
    row = addRow( p, row, "Fonte virtual: cabeçalho", comboRef, new JLabel(" = "), textRefVal );
    row = addRow( p, row, "Parear pelo cabeçalho", comboPair, new JLabel("(mesmo componente)"), null );
    row = addRow( p, row, "Normalização temporal", comboNorm, new JLabel(" janela [s]"), textRam );
    row = addRow( p, row, "Banda [Hz]  fmin", textFmin, new JLabel(" fmax"), textFmax );
    row = addRow( p, row, "", boxWhiten, null, null );
    row = addRow( p, row, "Lag máximo [s]", textLag, new JLabel(" janelas p/ empilhar [s]"), textWin );
    row = addRow( p, row, "", null, new JLabel(" sobreposição [%]"), textOverlap );
    row = addRow( p, row, "", boxSym, null, null );
    row = addRow( p, row, "", boxNormOut, null, null );
    row = addRow( p, row, "", boxOffset, null, null );
    row = addRow( p, row, "", boxDisp, null, null );
    row = addRow( p, row, "", boxAsym, null, null );
    JPanel asRow = new JPanel( new java.awt.FlowLayout( java.awt.FlowLayout.LEFT, 4, 0 ) );
    asRow.add( new JLabel( "janela de chegada: c mín" ) ); asRow.add( textAsCmin );
    asRow.add( new JLabel( "c máx [m/s]" ) ); asRow.add( textAsCmax );
    asRow.add( new JLabel( "folga [s]" ) ); asRow.add( textAsPad );
    row = addRow( p, row, "", asRow, null, null );

    int option = JOptionPane.showConfirmDialog( ctx.getParentFrame(), p, "Interferometria (VSG) - " + ctx.getTitle(),
                                                JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE );
    if( option != JOptionPane.OK_OPTION ) return;

    //--- Read parameters
    try {
      myParams.maxLag = Double.parseDouble( textLag.getText().trim() );
      myParams.fmin   = Double.parseDouble( textFmin.getText().trim() );
      myParams.fmax   = Double.parseDouble( textFmax.getText().trim() );
      myParams.ramWin = Double.parseDouble( textRam.getText().trim() );
      myParams.winLen = Double.parseDouble( textWin.getText().trim() );
      myParams.overlap = Double.parseDouble( textOverlap.getText().trim() );
      myAsymParams.cmin = Double.parseDouble( textAsCmin.getText().trim().replace( ',', '.' ) );
      myAsymParams.cmax = Double.parseDouble( textAsCmax.getText().trim().replace( ',', '.' ) );
      myAsymParams.pad  = Double.parseDouble( textAsPad.getText().trim().replace( ',', '.' ) );
    }
    catch( NumberFormatException e ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "Número inválido", "Interferometria", JOptionPane.ERROR_MESSAGE );
      return;
    }
    myParams.tempNorm  = comboNorm.getSelectedIndex();
    myParams.whiten    = boxWhiten.isSelected();
    myParams.symmetric = boxSym.isSelected();
    myParams.normalize = boxNormOut.isSelected();
    myRefHdr  = (String)comboRef.getSelectedItem();
    myRefVal  = textRefVal.getText().trim();
    myPairHdr = (String)comboPair.getSelectedItem();
    boolean setOffset = boxOffset.isSelected() && hasXY;
    myDisp = boxDisp.isSelected() && hasXY;
    myAsym = boxAsym.isSelected();
    if( myAsym && ( myAsymParams.cmin <= 0 || myAsymParams.cmax <= myAsymParams.cmin ) ) { error( ctx, "Assimetria: 0 < c mín < c máx" ); return; }

    if( myParams.maxLag <= 0 ) { error( ctx, "O lag máximo deve ser > 0" ); return; }
    if( myParams.overlap < 0 || myParams.overlap > 90 ) { error( ctx, "A sobreposição deve estar entre 0 e 90%" ); return; }
    if( myParams.fmax <= myParams.fmin ) { error( ctx, "fmax deve ser > fmin" ); return; }

    //--- Reference trace(s)
    int ihRef  = ctx.getHeaderIndex( myRefHdr );
    int ihPair = NONE.equals( myPairHdr ) ? -1 : ctx.getHeaderIndex( myPairHdr );
    int ntr = buf.numTraces();
    int[] refOf = new int[ntr];
    Arrays.fill( refOf, -1 );
    java.util.Map<String,Integer> refByPair = new java.util.LinkedHashMap<String,Integer>();
    for( int i = 0; i < ntr; i++ ) {
      if( !matches( buf.headerValues(i)[ihRef], myRefVal ) ) continue;
      String key = ( ihPair >= 0 ) ? buf.headerValues(i)[ihPair].toString() : "";
      if( !refByPair.containsKey( key ) ) refByPair.put( key, i );
    }
    if( refByPair.isEmpty() ) { error( ctx, "Nenhum traço exibido com " + myRefHdr + " = " + myRefVal ); return; }
    int numNoRef = 0;
    for( int i = 0; i < ntr; i++ ) {
      String key = ( ihPair >= 0 ) ? buf.headerValues(i)[ihPair].toString() : "";
      Integer r = refByPair.get( key );
      if( r != null ) refOf[i] = r; else numNoRef++;
    }

    //--- Compute
    float[][] samples = new float[ntr][];
    for( int i = 0; i < ntr; i++ ) samples[i] = buf.samples( i );
    Result res;
    csAssimetria.Result asym = null;
    if( myAsym ) {
      // asymmetry needs both sides, before summing and normalizing
      Params p2 = copy( myParams );
      p2.symmetric = false;
      p2.normalize = false;
      res = compute( samples, refOf, dt, p2 );
      double[] xs = null, ys = null;
      if( hasXY ) {
        int ix = ctx.getHeaderIndex( "rec_x" ), iy = ctx.getHeaderIndex( "rec_y" );
        xs = new double[ntr]; ys = new double[ntr];
        for( int i = 0; i < ntr; i++ ) { xs[i] = buf.headerValues(i)[ix].doubleValue(); ys[i] = buf.headerValues(i)[iy].doubleValue(); }
      }
      String[] gNames = refByPair.keySet().toArray( new String[0] );
      java.util.List<Integer> gRefs = new java.util.ArrayList<Integer>( refByPair.values() );
      int[] group = new int[ntr];
      for( int i = 0; i < ntr; i++ ) group[i] = refOf[i] >= 0 ? gRefs.indexOf( refOf[i] ) : -1;
      asym = csAssimetria.compute( res.corr, res.nl, dt, refOf, xs, ys, myAsymParams, gNames, group );
      for( int i = 0; i < ntr; i++ ) res.corr[i] = finish( res.corr[i], res.nl, myParams.symmetric, myParams.normalize );
    }
    else {
      res = compute( samples, refOf, dt, myParams );
    }

    //--- Output pane
    int nout = res.corr[0].length;
    int ihOff = ctx.getHeaderIndex("offset"), ihX = ctx.getHeaderIndex("rec_x"), ihY = ctx.getHeaderIndex("rec_y");
    csHeaderDef[] defsOut = defs;
    if( setOffset && ihOff < 0 ) {
      defsOut = Arrays.copyOf( defs, defs.length + 1 );
      defsOut[defs.length] = new csHeaderDef( "offset", "Distance to the virtual source [m] (Interferometria)", cseis.jni.csJNIDef.TYPE_DOUBLE );
      ihOff = defs.length;
    }
    int ihAsym = -1;
    if( asym != null ) {
      ihAsym = -1;
      for( int ih = 0; ih < defsOut.length; ih++ ) if( "asym".equals( defsOut[ih].name ) ) ihAsym = ih;
      if( ihAsym < 0 ) {
        defsOut = Arrays.copyOf( defsOut, defsOut.length + 1 );
        defsOut[defsOut.length - 1] = new csHeaderDef( "asym", "Causal/acausal asymmetry (E+ - E-)/(E+ + E-) (Interferometria)", cseis.jni.csJNIDef.TYPE_DOUBLE );
        ihAsym = defsOut.length - 1;
      }
    }
    csVirtualSeismicReader reader = new csVirtualSeismicReader( nout, defsOut.length, ctx.getSampleInt(), defsOut, ctx.getVerticalDomain() );
    csTraceBuffer out = reader.retrieveTraceBuffer();
    double[] dist = new double[ntr];
    Arrays.fill( dist, Double.NaN );
    for( int i = 0; i < ntr; i++ ) {
      csHeader[] hin = buf.headerValues( i );
      if( hasXY && refOf[i] >= 0 ) {
        csHeader[] hr = buf.headerValues( refOf[i] );
        dist[i] = Math.hypot( hin[ihX].doubleValue() - hr[ihX].doubleValue(), hin[ihY].doubleValue() - hr[ihY].doubleValue() );
      }
      csHeader[] hout = new csHeader[defsOut.length];
      for( int ih = 0; ih < defsOut.length; ih++ ) hout[ih] = ( ih < hin.length && ih < defs.length ) ? new csHeader( hin[ih] ) : new csHeader( 0 );
      if( setOffset && refOf[i] >= 0 ) hout[ihOff].setValue( dist[i] );
      if( ihAsym >= 0 ) hout[ihAsym].setValue( Double.isNaN( asym.asym[i] ) ? 0.0 : asym.asym[i] );
      out.addTrace( res.corr[i], hout );
    }
    ourCounter++;
    String title = String.format( Locale.US, "VSG%d %s (fonte %s=%s, %s%s%s, %.4g-%.4g Hz, %s)", ourCounter, ctx.getTitle(), myRefHdr, myRefVal,
        myParams.tempNorm == 0 ? "" : NORMS[myParams.tempNorm] + ", ", myParams.whiten ? "branqueado, " : "",
        myParams.symmetric ? "lados somados" : "dois lados", myParams.fmin, myParams.fmax,
        myParams.symmetric ? "lag 0 em 0 ms" : String.format( Locale.US, "lag 0 em %.0f ms", res.nl * dt * 1000.0 ) );
    ctx.openNewPane( reader, title );

    StringBuilder msg = new StringBuilder();
    msg.append( String.format( Locale.US, "%d correlações, %d janela(s) de %.4g s empilhada(s) (início a cada %.4g s, sobreposição %.0f%%), lag máximo %.4g s.%n",
        ntr, res.nWindows, res.winSamples * dt, res.hopSamples * dt, res.nWindows > 1 ? 100.0 * ( 1.0 - (double)res.hopSamples / res.winSamples ) : 0.0, res.nl * dt ) );
    msg.append( myParams.symmetric ? "Lados somados: eixo de tempo = lag 0 ... lag máximo.\n"
                                   : String.format( Locale.US, "Dois lados: o lag 0 está em %.0f ms no eixo de tempo (lags negativos acima).%n", res.nl * dt * 1000.0 ) );
    msg.append( "Fonte(s) virtual(is): " );
    for( java.util.Map.Entry<String,Integer> e : refByPair.entrySet() ) {
      msg.append( "traço " ).append( e.getValue()+1 ).append( ihPair >= 0 ? " (" + myPairHdr + " " + e.getKey() + ")" : "" ).append( "  " );
    }
    if( numNoRef > 0 ) msg.append( String.format( "%n%d traços sem referência com o mesmo %s: saída zerada.", numNoRef, myPairHdr ) );
    if( setOffset ) msg.append( "\nCabeçalho 'offset' = distância à fonte virtual (use Sort by header em offset)." );
    if( asym != null ) msg.append( "\nCabeçalho 'asym' = assimetria causal × acausal (0 na fonte virtual)." );
    JOptionPane.showMessageDialog( ctx.getParentFrame(), msg.toString(), "Interferometria (VSG)", JOptionPane.INFORMATION_MESSAGE );

    //--- asymmetry window
    if( asym != null ) {
      int ihTrc = ctx.getHeaderIndex( "trcno" );
      double[] trc = new double[ntr];
      for( int i = 0; i < ntr; i++ ) trc[i] = ihTrc >= 0 ? buf.headerValues(i)[ihTrc].doubleValue() : i + 1;
      csAssimetria.Frame af = new csAssimetria.Frame( asym, myAsymParams, "VSG" + ourCounter + " " + ctx.getTitle(), trc );
      af.setLocationRelativeTo( ctx.getParentFrame() );
      af.setVisible( true );
    }

    //--- f-c dispersion image of the VSG (one image per virtual source / pairing value)
    if( myDisp ) {
      for( java.util.Map.Entry<String,Integer> e : refByPair.entrySet() ) {
        java.util.List<Integer> idx = new java.util.ArrayList<Integer>();
        for( int i = 0; i < ntr; i++ ) if( refOf[i] == e.getValue() ) idx.add( i );
        float[][] tr = new float[idx.size()][];
        double[] off = new double[idx.size()];
        for( int k = 0; k < tr.length; k++ ) { tr[k] = res.corr[idx.get( k )]; off[k] = dist[idx.get( k )]; }
        csDispersionFC.Params dp = new csDispersionFC.Params();
        dp.oneSided = myParams.symmetric;
        dp.lag0 = myParams.symmetric ? 0 : res.nl;
        dp.side = csDispersionFC.SIDE_SUM;
        dp.fmin = Math.max( myParams.fmin, 0.0 );
        dp.fmax = Math.min( myParams.fmax, 0.5 / dt );
        dp.dmin = 1.0;                       // the virtual source itself (offset 0) carries no velocity information
        String t = "VSG" + ourCounter + " " + ctx.getTitle() + " (fonte " + myRefHdr + "=" + myRefVal
                   + ( ihPair >= 0 ? ", " + myPairHdr + " " + e.getKey() : "" ) + ")";
        csPluginDispersao.show( ctx, tr, off, dt, dp, t, "distância à fonte virtual" );
      }
    }
  }

  private static Params copy( Params a ) {
    Params b = new Params();
    b.maxLag = a.maxLag; b.fmin = a.fmin; b.fmax = a.fmax; b.tempNorm = a.tempNorm; b.ramWin = a.ramWin;
    b.whiten = a.whiten; b.winLen = a.winLen; b.overlap = a.overlap; b.symmetric = a.symmetric; b.normalize = a.normalize;
    return b;
  }
  /** Two-sided correlation (lag 0 at nl) -> output trace: optional sum of both sides, optional max normalization */
  static float[] finish( float[] two, int nl, boolean symmetric, boolean normalize ) {
    float[] o;
    if( symmetric ) {
      o = new float[nl + 1];
      for( int j = 0; j <= nl; j++ ) o[j] = two[nl + j] + two[nl - j];
    }
    else o = two.clone();
    if( normalize ) {
      float m = 0f;
      for( float v : o ) m = Math.max( m, Math.abs( v ) );
      if( m > 0f ) for( int j = 0; j < o.length; j++ ) o[j] /= m;
    }
    return o;
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
  private static boolean matches( csHeader h, String text ) {
    try {
      return Math.abs( h.doubleValue() - Double.parseDouble( text ) ) < 1.0e-6 * Math.max( 1.0, Math.abs( h.doubleValue() ) );
    }
    catch( NumberFormatException e ) {
      return h.toString().trim().equals( text );
    }
  }
  private static void error( csIPluginContext ctx, String text ) {
    JOptionPane.showMessageDialog( ctx.getParentFrame(), text, "Interferometria", JOptionPane.ERROR_MESSAGE );
  }
  private static String fmt( double v ) {
    if( v == Math.rint( v ) ) return String.valueOf( (long)v );
    return String.format( Locale.US, "%.4g", v );
  }

  //====================================================================
  // Computation (no Swing)
  //====================================================================
  public static final class Result {
    public float[][] corr;    // [trace][sample]
    public int nl;            // max lag in samples
    public int nWindows;
    public int winSamples;
    public int hopSamples;
  }

  /**
   * @param samples traces [ntr][ns]
   * @param refOf   index of the reference trace for each trace (-1: none, output zero)
   * @param dt      sample interval [s]
   */
  public static Result compute( float[][] samples, int[] refOf, double dt, Params p ) {
    int ntr = samples.length;
    int ns = samples[0].length;
    Result r = new Result();
    int nw = ( p.winLen > 0 ) ? Math.min( ns, (int)Math.round( p.winLen / dt ) ) : ns;
    nw = Math.max( nw, 4 );
    int nl = Math.min( (int)Math.round( p.maxLag / dt ), nw - 1 );
    // Windows start every 'hop' samples (overlap in %); a window must fit completely in the trace
    int hop = Math.max( 1, (int)Math.round( nw * ( 1.0 - Math.max( 0.0, Math.min( 90.0, p.overlap ) ) / 100.0 ) ) );
    if( nw >= ns ) hop = nw;
    int nwin = Math.max( 1, 1 + ( ns - nw ) / hop );
    int nfft = 1;
    while( nfft < nw + nl ) nfft *= 2;
    r.nl = nl; r.nWindows = nwin; r.winSamples = nw; r.hopSamples = hop;

    double df = 1.0 / ( nfft * dt );
    double fmax = ( p.fmax > 0 ) ? Math.min( p.fmax, 0.5/dt ) : 0.5/dt;
    double[] mask = bandMask( nfft, df, p.fmin, fmax );
    double[] taper = tukey( nw, 0.05 );
    int nram = Math.max( 1, (int)Math.round( p.ramWin / dt ) );

    // Only the spectra of the reference traces are kept (all windows); every other trace is
    // pre-processed window by window and discarded, so memory does not grow with the number of traces
    boolean[] isRef = new boolean[ntr];
    for( int i = 0; i < ntr; i++ ) if( refOf[i] >= 0 ) isRef[refOf[i]] = true;
    double[][][] refRe = new double[ntr][][], refIm = new double[ntr][][];
    for( int i = 0; i < ntr; i++ ) {
      if( !isRef[i] ) continue;
      refRe[i] = new double[nwin][];
      refIm[i] = new double[nwin][];
      for( int w = 0; w < nwin; w++ ) {
        double[] re = new double[nfft], im = new double[nfft];
        preprocess( samples[i], w * hop, nw, taper, mask, nram, p, re, im );
        refRe[i][w] = re;
        refIm[i][w] = im;
      }
    }

    int nout = p.symmetric ? nl + 1 : 2 * nl + 1;
    r.corr = new float[ntr][nout];
    double[] re = new double[nfft], im = new double[nfft], acc = new double[2*nl+1];
    double[] xr = new double[nfft], xi = new double[nfft];
    for( int i = 0; i < ntr; i++ ) {
      int k = refOf[i];
      if( k < 0 ) continue;
      Arrays.fill( acc, 0.0 );
      for( int w = 0; w < nwin; w++ ) {
        double[] rr = refRe[k][w], ri = refIm[k][w];
        if( isRef[i] ) {
          System.arraycopy( refRe[i][w], 0, xr, 0, nfft );
          System.arraycopy( refIm[i][w], 0, xi, 0, nfft );
        }
        else {
          preprocess( samples[i], w * hop, nw, taper, mask, nram, p, xr, xi );
        }
        for( int f = 0; f < nfft; f++ ) {          // conj(R) * X
          re[f] = rr[f] * xr[f] + ri[f] * xi[f];
          im[f] = rr[f] * xi[f] - ri[f] * xr[f];
        }
        csClockXcorr.ifft( re, im );
        for( int j = -nl; j <= nl; j++ ) acc[j + nl] += re[ ( j + nfft ) % nfft ];
      }
      float[] o = r.corr[i];
      if( p.symmetric ) {
        for( int j = 0; j <= nl; j++ ) o[j] = (float)( acc[nl + j] + acc[nl - j] );
      }
      else {
        for( int j = 0; j < 2*nl+1; j++ ) o[j] = (float)acc[j];
      }
      if( p.normalize ) {
        float m = 0f;
        for( float v : o ) m = Math.max( m, Math.abs( v ) );
        if( m > 0f ) for( int j = 0; j < o.length; j++ ) o[j] /= m;
      }
    }
    return r;
  }

  static double[] bandMask( int nfft, double df, double fmin, double fmax ) {
    double[] m = new double[nfft];
    double w = Math.max( 2*df, 0.1 * ( fmax - fmin ) );   // cosine transition width
    for( int f = 0; f <= nfft/2; f++ ) {
      double fr = f * df, v;
      if( fr < fmin - w || fr > fmax + w ) v = 0.0;
      else if( fr < fmin ) v = ( fmin <= 0 ) ? 1.0 : 0.5 - 0.5 * Math.cos( Math.PI * ( fr - ( fmin - w ) ) / w );
      else if( fr > fmax ) v = 0.5 + 0.5 * Math.cos( Math.PI * ( fr - fmax ) / w );
      else v = 1.0;
      if( fmin <= 0 && fr <= fmin ) v = ( f == 0 ) ? 0.0 : v;   // no DC
      m[f] = v;
      if( f > 0 && f < nfft/2 ) m[nfft - f] = v;
    }
    return m;
  }
  static double[] tukey( int n, double frac ) {
    double[] w = new double[n];
    int ne = Math.max( 1, (int)( frac * n ) );
    for( int i = 0; i < n; i++ ) {
      double v = 1.0;
      if( i < ne ) v = 0.5 - 0.5 * Math.cos( Math.PI * i / ne );
      else if( i >= n - ne ) v = 0.5 - 0.5 * Math.cos( Math.PI * ( n - 1 - i ) / ne );
      w[i] = v;
    }
    return w;
  }
  /** Window [i0, i0+nw) of trace s -> pre-processed spectrum (length nfft = re.length) */
  static void preprocess( float[] s, int i0, int nw, double[] taper, double[] mask, int nram, Params p, double[] re, double[] im ) {
    int nfft = re.length;
    Arrays.fill( re, 0.0 );
    Arrays.fill( im, 0.0 );
    double mean = 0.0;
    for( int i = 0; i < nw; i++ ) { re[i] = s[i0 + i]; mean += re[i]; }
    mean /= nw;
    for( int i = 0; i < nw; i++ ) re[i] = ( re[i] - mean ) * taper[i];
    csClockXcorr.fft( re, im );
    for( int f = 0; f < nfft; f++ ) { re[f] *= mask[f]; im[f] *= mask[f]; }
    if( p.tempNorm != 0 ) {
      csClockXcorr.ifft( re, im );
      if( p.tempNorm == 1 ) {                         // one-bit
        for( int i = 0; i < nfft; i++ ) re[i] = ( i < nw ) ? Math.signum( re[i] ) * taper[i] : 0.0;
      }
      else {                                          // running absolute mean
        double[] a = new double[nw + 1];
        for( int i = 0; i < nw; i++ ) a[i + 1] = a[i] + Math.abs( re[i] );
        int h = nram / 2;
        double[] o = new double[nw];
        for( int i = 0; i < nw; i++ ) {
          int j0 = Math.max( 0, i - h ), j1 = Math.min( nw, i + h + 1 );
          double m = ( a[j1] - a[j0] ) / ( j1 - j0 );
          o[i] = ( m > 0 ) ? re[i] / m : 0.0;
        }
        for( int i = 0; i < nfft; i++ ) re[i] = ( i < nw ) ? o[i] * taper[i] : 0.0;
      }
      Arrays.fill( im, 0.0 );
      csClockXcorr.fft( re, im );
      for( int f = 0; f < nfft; f++ ) { re[f] *= mask[f]; im[f] *= mask[f]; }
    }
    if( p.whiten ) {
      // divide by the amplitude spectrum smoothed over a few bins, inside the band
      int h = nfft / 2 + 1;
      double[] amp = new double[h];
      for( int f = 0; f < h; f++ ) amp[f] = Math.hypot( re[f], im[f] );
      int half = Math.max( 1, nfft / 400 );
      double[] c = new double[h + 1];
      for( int f = 0; f < h; f++ ) c[f + 1] = c[f] + amp[f];
      double big = 0.0;
      double[] sm = new double[h];
      for( int f = 0; f < h; f++ ) {
        int f0 = Math.max( 0, f - half ), f1 = Math.min( h, f + half + 1 );
        sm[f] = ( c[f1] - c[f0] ) / ( f1 - f0 );
        big = Math.max( big, sm[f] );
      }
      double eps = 1.0e-6 * big;
      for( int f = 0; f < h; f++ ) {
        double g = mask[f] / ( sm[f] + eps );
        re[f] *= g; im[f] *= g;
        if( f > 0 && f < nfft / 2 ) { re[nfft - f] *= g; im[nfft - f] *= g; }
      }
    }
  }
}
