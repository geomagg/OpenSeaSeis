/* SeaView plugin: clock QC by noise cross-correlation. M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import cseis.general.csUnits;
import cseis.jni.csVirtualSeismicReader;
import cseis.seis.csHeader;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicTraceBuffer;
import cseis.seis.csTraceBuffer;
import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.LinkedHashSet;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ProgressMonitor;
import javax.swing.SwingWorker;

/**
 * Plugin: clock QC of the nodes of the active pane (typically one receiver line,
 * LxxxxNyyyy-Nzzzz.cseis, one trace per node/sensor) by ambient-noise cross-correlation.
 * <p>
 * 1) one trace per node (chosen sensor), 2) alignment of all traces to the common absolute time
 * window, 3) pre-processing + cross-correlation of each node with its K next neighbours, stacked,
 * 4) clock error from the time symmetry of each correlation (positive x negative lags), solved per
 * node by least squares, also per period (clock drift). See csClockXcorr for the method.
 * <p>
 * Output: a new pane with the stacked correlations (one trace per pair, lag 0 at the centre of the
 * trace) and a report window (per node: clock error, drift; per pair: shift, quality, A+/A-),
 * which can be saved as CSV.
 */
public class csPluginClockXcorr implements csISeaViewPlugin {
  private static int ourCounter = 0;
  private final csClockXcorr.Params myParams = new csClockXcorr.Params();
  private String myNodeHdr = null;
  private String mySensorHdr = null;
  private String mySensorValue = null;
  private String myRefNode = "";

  @Override
  public String getName() {
    return "Clock: correlação cruzada dos nodes...";
  }
  @Override
  public String getDescription() {
    return "Alinha os traços no tempo absoluto, correlaciona cada node com os vizinhos e mede o erro de clock pela assimetria (lags positivos x negativos) da correlação";
  }

  @Override
  public void run( csIPluginContext ctx ) {
    final String title = "Clock (correlação cruzada)";
    csISeismicTraceBuffer buf = ctx.getTraceBuffer();
    if( buf == null || buf.numTraces() < 2 ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "São necessários pelo menos 2 traços no painel ativo", title, JOptionPane.WARNING_MESSAGE );
      return;
    }
    if( ctx.getVerticalDomain() != csUnits.DOMAIN_TIME ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "O painel ativo não está no domínio do tempo", title, JOptionPane.WARNING_MESSAGE );
      return;
    }
    final csHeaderDef[] defs = ctx.getHeaderDef();

    //--- Node / sensor headers: detected automatically (no question). The "Headers..." button of
    //    the dialog shows every header of the 1st trace and lets you change them if needed.
    final int[] hdrSel = { detectHeader( defs, myNodeHdr, NODE_HDRS ), detectHeader( defs, mySensorHdr, SENSOR_HDRS ) };
    if( hdrSel[0] < 0 ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "Nenhum header de node encontrado ("
          + String.join( ", ", NODE_HDRS ) + ").\nHeaders do 1º traço:\n" + headerList( defs, buf.headerValues( 0 ), 40 ),
          title, JOptionPane.ERROR_MESSAGE );
      return;
    }
    JLabel infoHdr = new JLabel();
    JComboBox<String> comboSensorVal = new JComboBox<>();
    Runnable refreshHdr = () -> {
      List<String> vals = sensorValues( buf, hdrSel[0], hdrSel[1] );
      comboSensorVal.removeAllItems();
      for( String v : vals ) comboSensorVal.addItem( v );
      if( mySensorValue != null && vals.contains( mySensorValue ) ) comboSensorVal.setSelectedItem( mySensorValue );
      comboSensorVal.setEnabled( vals.size() > 1 );
      String sensorDesc = ( hdrSel[1] >= 0 ) ? "header <b>" + defs[hdrSel[1]].name + "</b>"
                                             : "sem header de sensor: <b>ordem do traço dentro do node</b>";
      infoHdr.setText( "<html>Node: header <b>" + defs[hdrSel[0]].name + "</b> (" + countNodes( buf, hdrSel[0] )
          + " nodes) &nbsp; | &nbsp; Sensor: " + sensorDesc + " (" + vals.size() + " valor(es))</html>" );
    };
    refreshHdr.run();
    JButton btnHdr = new JButton( "Headers..." );
    btnHdr.setToolTipText( "Ver todos os headers do 1º traço e, se precisar, trocar os headers de node/sensor" );
    btnHdr.addActionListener( e -> {
      if( chooseHeaders( ctx, defs, buf, hdrSel ) ) refreshHdr.run();
    } );

    //--- Parameter dialog
    JTextField textFmin = num( myParams.fmin ), textFmax = num( myParams.fmax );
    JTextField textWin = num( myParams.winSec ), textLag = num( myParams.maxLagSec ), textMinLag = num( myParams.minLagSec );
    JTextField textPer = num( myParams.periodMin );
    JTextField textK = new JTextField( "" + myParams.neighbors );
    JTextField textQ = num( myParams.minQuality );
    JTextField textRef = new JTextField( myRefNode );
    textRef.setToolTipText( "Node de referência (erro = 0). Vazio: média dos erros = 0" );
    JCheckBox boxOneBit = new JCheckBox( "Normalização one-bit", myParams.oneBit );
    JCheckBox boxWhiten = new JCheckBox( "Branqueamento espectral", myParams.whiten );

    JPanel panel = new JPanel( new GridLayout( 0, 2, 6, 4 ) );
    panel.add( new JLabel( "Sensor a usar:" ) );                      panel.add( comboSensorVal );
    panel.add( new JLabel( "Banda: freq. mínima [Hz]:" ) );           panel.add( textFmin );
    panel.add( new JLabel( "Banda: freq. máxima [Hz]:" ) );           panel.add( textFmax );
    panel.add( new JLabel( "Janela de correlação [s]:" ) );           panel.add( textWin );
    panel.add( new JLabel( "Lag máximo [s]:" ) );                     panel.add( textLag );
    panel.add( new JLabel( "Ignorar |lag| menor que [s]:" ) );        panel.add( textMinLag );
    panel.add( new JLabel( "Período p/ deriva [min] (0 = só total):" ) ); panel.add( textPer );
    panel.add( new JLabel( "Pares: vizinhos por node (K):" ) );       panel.add( textK );
    panel.add( new JLabel( "Qualidade mínima do par (0-1):" ) );      panel.add( textQ );
    panel.add( new JLabel( "Node de referência (vazio = média):" ) ); panel.add( textRef );
    panel.add( boxOneBit );                                           panel.add( boxWhiten );
    JPanel top = new JPanel( new BorderLayout( 6, 6 ) );
    top.add( new JLabel( "<html>Antes da correlação todos os traços são alinhados no tempo absoluto:<br>"
        + "começam no início mais tardio e terminam no fim mais cedo (headers time_samp1<br>"
        + "ou time_year/day/hour/min/sec).</html>" ), BorderLayout.NORTH );
    top.add( infoHdr, BorderLayout.CENTER );
    top.add( btnHdr, BorderLayout.EAST );
    JPanel outer = new JPanel( new BorderLayout( 4, 8 ) );
    outer.add( top, BorderLayout.NORTH );
    outer.add( panel, BorderLayout.CENTER );
    int option = JOptionPane.showConfirmDialog( ctx.getParentFrame(), outer, title + " - " + ctx.getTitle(),
        JOptionPane.OK_CANCEL_OPTION, JOptionPane.PLAIN_MESSAGE );
    if( option != JOptionPane.OK_OPTION ) return;
    try {
      myParams.fmin = dbl( textFmin );
      myParams.fmax = dbl( textFmax );
      myParams.winSec = dbl( textWin );
      myParams.maxLagSec = dbl( textLag );
      myParams.minLagSec = dbl( textMinLag );
      myParams.periodMin = dbl( textPer );
      myParams.neighbors = Integer.parseInt( textK.getText().trim() );
      myParams.minQuality = dbl( textQ );
      myRefNode = textRef.getText().trim();
      myParams.refNode = myRefNode.isEmpty() ? Integer.MIN_VALUE : Integer.parseInt( myRefNode );
    }
    catch( NumberFormatException e ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "Número inválido: " + e.getMessage(), title, JOptionPane.ERROR_MESSAGE );
      return;
    }
    myParams.oneBit = boxOneBit.isSelected();
    myParams.whiten = boxWhiten.isSelected();
    final int iNode = hdrSel[0], iSensor = hdrSel[1];
    myNodeHdr = defs[iNode].name;
    mySensorHdr = ( iSensor >= 0 ) ? defs[iSensor].name : null;
    final String sensorUsed = (String)comboSensorVal.getSelectedItem();
    mySensorValue = sensorUsed;
    double dt = ctx.getSampleInt() / 1000.0;
    if( !( myParams.fmin >= 0.0 && myParams.fmax > myParams.fmin && myParams.fmax < 0.5 / dt ) ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), String.format( Locale.US,
          "Banda inválida (precisa 0 <= fmin < fmax < Nyquist = %.4g Hz)", 0.5 / dt ), title, JOptionPane.ERROR_MESSAGE );
      return;
    }

    //--- One trace per node (the chosen sensor)
    Map<Integer,Integer> traceOfNode = new LinkedHashMap<>();
    Map<Integer,Integer> countInNode = new LinkedHashMap<>();
    int nDup = 0, nNoTime = 0;
    for( int i = 0; i < buf.numTraces(); i++ ) {
      csHeader[] h = buf.headerValues( i );
      int node = h[iNode].intValue();
      int order = countInNode.merge( node, 1, Integer::sum );
      String sv = ( iSensor >= 0 ) ? hdrText( h[iSensor] ) : String.valueOf( order );
      if( sensorUsed != null && !sameValue( sv, sensorUsed ) ) continue;
      if( Double.isNaN( csClockXcorr.startTimeSec( h, defs ) ) ) { nNoTime++; continue; }
      if( traceOfNode.containsKey( node ) ) { nDup++; continue; }
      traceOfNode.put( node, i );
    }
    if( traceOfNode.size() < 2 ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), String.format( Locale.US,
          "Menos de 2 nodes utilizáveis com o sensor '%s'.%nTraços sem tempo no header: %d%n%nHeaders do 1º traço:%n%s",
          sensorUsed, nNoTime, headerList( defs, buf.headerValues( 0 ), 40 ) ), title, JOptionPane.ERROR_MESSAGE );
      return;
    }
    final int n = traceOfNode.size();
    final int[] nodes = new int[n];
    final float[][] samples = new float[n][];
    final double[] t0 = new double[n];
    int k = 0;
    for( Map.Entry<Integer,Integer> en : traceOfNode.entrySet() ) {
      nodes[k] = en.getKey();
      int itr = en.getValue();
      samples[k] = buf.samples( itr );
      t0[k] = csClockXcorr.startTimeSec( buf.headerValues( itr ), defs );
      k++;
    }
    final csClockXcorr.Params params = copy( myParams );
    final String paneTitle = ctx.getTitle();
    final String sensorTxt = ( iSensor >= 0 ) ? mySensorHdr + "=" + sensorUsed : "sensor (ordem no node)=" + sensorUsed;
    final int dupCount = nDup, noTimeCount = nNoTime;

    //--- Compute in background
    final ProgressMonitor monitor = new ProgressMonitor( ctx.getParentFrame(), title + " - " + paneTitle, "Preparando...", 0, 1000 );
    monitor.setMillisToDecideToPopup( 200 );
    monitor.setMillisToPopup( 200 );
    SwingWorker<csClockXcorr.Result,Void> worker = new SwingWorker<csClockXcorr.Result,Void>() {
      private String myError = null;
      @Override
      protected csClockXcorr.Result doInBackground() {
        try {
          return csClockXcorr.run( nodes, samples, t0, dt, params, ( done, total, msg ) -> {
            final int v = (int)( 1000L * done / Math.max( 1, total ) );
            javax.swing.SwingUtilities.invokeLater( () -> { monitor.setProgress( v ); monitor.setNote( msg ); } );
            return !monitor.isCanceled();
          } );
        }
        catch( Throwable t ) {
          myError = t.getMessage() != null ? t.getMessage() : t.toString();
          return null;
        }
      }
      @Override
      protected void done() {
        monitor.close();
        csClockXcorr.Result res;
        try { res = get(); } catch( Exception e ) { res = null; myError = e.toString(); }
        if( res == null ) {
          if( myError != null ) JOptionPane.showMessageDialog( ctx.getParentFrame(), myError, title, JOptionPane.ERROR_MESSAGE );
          return;
        }
        ourCounter++;
        openCorrelationPane( ctx, res, "XCLK" + ourCounter + " " + paneTitle );
        showReport( ctx, res, params, paneTitle, sensorTxt, dupCount, noTimeCount );
      }
    };
    worker.execute();
  }

  //--------------------------------------------------------------------
  /** New pane: one trace per pair, 2*nl+1 samples, lag 0 at sample nl */
  private static void openCorrelationPane( csIPluginContext ctx, csClockXcorr.Result r, String title ) {
    csHeaderDef[] defs = {
      new csHeaderDef( "trcno", "Trace number", csHeaderDef.TYPE_INT ),
      new csHeaderDef( "node_a", "First node of the pair", csHeaderDef.TYPE_INT ),
      new csHeaderDef( "node_b", "Second node of the pair", csHeaderDef.TYPE_INT ),
      new csHeaderDef( "rcv", "Second node of the pair (rcv)", csHeaderDef.TYPE_INT ),
      new csHeaderDef( "dnode", "node_b - node_a", csHeaderDef.TYPE_INT ),
      new csHeaderDef( "shift_ms", "Clock shift e_b - e_a [ms]", csHeaderDef.TYPE_FLOAT ),
      new csHeaderDef( "quality", "Symmetry correlation (0-1)", csHeaderDef.TYPE_FLOAT ),
      new csHeaderDef( "amp_ratio", "max|C| lags>0 / max|C| lags<0", csHeaderDef.TYPE_FLOAT ),
      new csHeaderDef( "lag0_ms", "Time of lag 0 in this trace [ms]", csHeaderDef.TYPE_FLOAT ),
    };
    int len = 2 * r.nl + 1;
    float dtMs = (float)( r.dt * 1000.0 );
    csVirtualSeismicReader reader = new csVirtualSeismicReader( len, defs.length, dtMs, defs, csUnits.DOMAIN_TIME );
    csTraceBuffer out = reader.retrieveTraceBuffer();
    for( int p = 0; p < r.pairA.length; p++ ) {
      float[] s = r.cc[p].clone();
      float mx = 0.0f;
      for( float v : s ) mx = Math.max( mx, Math.abs( v ) );
      if( mx > 0.0f ) for( int i = 0; i < len; i++ ) s[i] /= mx;
      int a = r.nodes[r.pairA[p]], b = r.nodes[r.pairB[p]];
      csHeader[] h = new csHeader[defs.length];
      for( int i = 0; i < h.length; i++ ) h[i] = new csHeader();
      h[0].setValue( p + 1 );
      h[1].setValue( a );
      h[2].setValue( b );
      h[3].setValue( b );
      h[4].setValue( b - a );
      h[5].setValue( (float)( 1000.0 * r.pairShift[p] ) );
      h[6].setValue( (float)r.pairQuality[p] );
      h[7].setValue( (float)( r.pairAmpNeg[p] > 0.0 ? r.pairAmpPos[p] / r.pairAmpNeg[p] : 0.0 ) );
      h[8].setValue( (float)( r.nl * r.dt * 1000.0 ) );
      out.addTrace( s, h );
    }
    ctx.openNewPane( reader, String.format( Locale.US, "%s (lag 0 = %.0f ms)", title, r.nl * r.dt * 1000.0 ) );
  }

  //--------------------------------------------------------------------
  private static void showReport( csIPluginContext ctx, csClockXcorr.Result r, csClockXcorr.Params p,
                                  String paneTitle, String sensorTxt, int nDup, int nNoTime ) {
    String text = report( r, p, paneTitle, sensorTxt, nDup, nNoTime );
    JTextArea area = new JTextArea( text, 34, 100 );
    area.setEditable( false );
    area.setFont( new Font( Font.MONOSPACED, Font.PLAIN, 12 ) );
    area.setCaretPosition( 0 );
    JDialog dialog = new JDialog( ctx.getParentFrame(), "Clock dos nodes - " + paneTitle, false );
    JButton save = new JButton( "Salvar CSV..." );
    save.addActionListener( e -> saveCsv( dialog, r ) );
    JButton close = new JButton( "Fechar" );
    close.addActionListener( e -> dialog.dispose() );
    JPanel buttons = new JPanel( new FlowLayout( FlowLayout.RIGHT ) );
    buttons.add( save );
    buttons.add( close );
    dialog.getContentPane().add( new JScrollPane( area ), BorderLayout.CENTER );
    dialog.getContentPane().add( buttons, BorderLayout.SOUTH );
    dialog.pack();
    dialog.setLocationRelativeTo( ctx.getParentFrame() );
    dialog.setVisible( true );
  }

  static String report( csClockXcorr.Result r, csClockXcorr.Params p, String paneTitle, String sensorTxt,
                        int nDup, int nNoTime ) {
    StringBuilder sb = new StringBuilder();
    csClockXcorr.Alignment al = r.align;
    sb.append( String.format( Locale.US, "Painel: %s   |   %s   |   %d nodes%n", paneTitle, sensorTxt, r.nodes.length ) );
    if( nDup > 0 ) sb.append( String.format( Locale.US, "  (%d traço(s) repetidos do mesmo node ignorados)%n", nDup ) );
    if( nNoTime > 0 ) sb.append( String.format( Locale.US, "  (%d traço(s) sem tempo no header ignorados)%n", nNoTime ) );
    sb.append( String.format( Locale.US, "%nALINHAMENTO (janela comum)%n" ) );
    sb.append( String.format( Locale.US, "  início: %s   fim: %s   duração: %.3f s (%d amostras, dt %.4g ms)%n",
        csProcessingAlignTime.fmt( al.tStart ), csProcessingAlignTime.fmt( al.tStart + ( al.nCommon - 1 ) * r.dt ),
        ( al.nCommon - 1 ) * r.dt, al.nCommon, r.dt * 1000.0 ) );
    sb.append( String.format( Locale.US, "  banda %.3g-%.3g Hz, janela %.4g s (%d janelas), lag máx %.4g s, |lag| >= %.4g s, K=%d, one-bit %s, branqueamento %s%n",
        p.fmin, p.fmax, p.winSec, r.nWindows, r.nl * r.dt, p.minLagSec, p.neighbors, p.oneBit ? "sim" : "não", p.whiten ? "sim" : "não" ) );
    sb.append( String.format( Locale.US, "  referência: %s%n", p.refNode == Integer.MIN_VALUE ? "média dos erros = 0" : "node " + p.refNode + " = 0" ) );

    sb.append( String.format( Locale.US, "%nPOR NODE  (erro > 0: relógio ATRASADO -- o evento aparece mais cedo no traço)%n" ) );
    sb.append( String.format( Locale.US, "  %8s %14s %12s %14s %14s %8s%n", "node", "inicio orig.[s]", "erro [ms]",
        "erro-tend.[ms]", "deriva[ms/dia]", "pares" ) );
    int nSusp = 0;
    for( int i = 0; i < r.nodes.length; i++ ) {
      if( r.nodeSuspect[i] ) nSusp++;
      sb.append( String.format( Locale.US, "  %8d %+14.3f %12s %14s %14s %8d%s%n", r.nodes[i], al.t0[i] - al.tStart,
          fmtMs( r.nodeErr[i] ), fmtMs( r.nodeResid[i] ), fmtMs( r.nodeDrift[i] ), r.nodePairsUsed[i],
          r.nodeSuspect[i] ? "   <== SUSPEITO" : "" ) );
    }
    sb.append( String.format( Locale.US, "  erro-tend. = erro menos a tendência linear ao longo da linha (ruído direcional gera%n"
        + "  uma tendência suave; um relógio com problema sai fora dela). %s%n",
        nSusp == 0 ? "Nenhum node suspeito." : nSusp + " node(s) suspeito(s)." ) );
    if( r.nPeriods > 1 ) {
      sb.append( String.format( Locale.US, "%nERRO POR PERÍODO [ms]  (centro do período, h após o início comum)%n  %8s", "node" ) );
      for( int q = 0; q < r.nPeriods; q++ ) sb.append( String.format( Locale.US, " %9.2fh", r.periodTime[q] / 3600.0 ) );
      sb.append( String.format( "%n" ) );
      for( int i = 0; i < r.nodes.length; i++ ) {
        sb.append( String.format( Locale.US, "  %8d", r.nodes[i] ) );
        for( int q = 0; q < r.nPeriods; q++ ) sb.append( String.format( Locale.US, " %10s", fmtMs( r.nodeErrPeriod[i][q] ) ) );
        sb.append( String.format( "%n" ) );
      }
    }
    sb.append( String.format( Locale.US, "%nPOR PAR  (deslocamento = e_b - e_a; qualidade = correlação lags>0 x lags<0; A+/A- = assimetria de amplitude)%n" ) );
    sb.append( String.format( Locale.US, "  %8s %8s %16s %10s %8s%n", "node_a", "node_b", "desloc. [ms]", "qualidade", "A+/A-" ) );
    for( int k = 0; k < r.pairA.length; k++ ) {
      sb.append( String.format( Locale.US, "  %8d %8d %16s %10.2f %8.2f%s%n", r.nodes[r.pairA[k]], r.nodes[r.pairB[k]],
          fmtMs( r.pairShift[k] ), r.pairQuality[k],
          r.pairAmpNeg[k] > 0.0 ? r.pairAmpPos[k] / r.pairAmpNeg[k] : 0.0,
          r.pairQuality[k] < p.minQuality ? "  (não usado)" : "" ) );
    }
    sb.append( String.format( Locale.US, "%nDeriva: inclinação do erro por período (só com períodos cobrindo >= %.0f min).%n",
        csClockXcorr.MIN_DRIFT_SPAN_SEC / 60.0 ) );
    sb.append( String.format( Locale.US, "Obs.: o erro só é confiável se for menor que o tempo de trânsito entre os nodes do par;%n"
        + "pares vizinhos muito próximos podem errar com erros grandes -- aumente K. Com 'média = 0', um node%n"
        + "com deriva puxa os outros em 1/n: informe um node de referência bom.%n" ) );
    return sb.toString();
  }
  private static String fmtMs( double sec ) {
    return Double.isNaN( sec ) ? "-" : String.format( Locale.US, "%+.2f", 1000.0 * sec );
  }

  private static void saveCsv( JDialog parent, csClockXcorr.Result r ) {
    JFileChooser fc = new JFileChooser();
    fc.setSelectedFile( new File( "clock_nodes.csv" ) );
    if( fc.showSaveDialog( parent ) != JFileChooser.APPROVE_OPTION ) return;
    File f = fc.getSelectedFile();
    try( PrintWriter w = new PrintWriter( new FileWriter( f ) ) ) {
      w.print( "node,start_offset_s,clock_error_ms,error_minus_trend_ms,suspect,drift_ms_per_day,pairs_used" );
      for( int q = 0; q < r.nPeriods; q++ ) w.printf( Locale.US, ",err_ms_%.3fh", r.periodTime[q] / 3600.0 );
      w.println();
      for( int i = 0; i < r.nodes.length; i++ ) {
        w.printf( Locale.US, "%d,%.6f,%s,%s,%d,%s,%d", r.nodes[i], r.align.t0[i] - r.align.tStart,
            csv( r.nodeErr[i] ), csv( r.nodeResid[i] ), r.nodeSuspect[i] ? 1 : 0, csv( r.nodeDrift[i] ), r.nodePairsUsed[i] );
        for( int q = 0; q < r.nPeriods; q++ ) w.print( "," + csv( r.nodeErrPeriod[i][q] ) );
        w.println();
      }
    }
    catch( IOException e ) {
      JOptionPane.showMessageDialog( parent, "Erro ao salvar: " + e.getMessage(), "CSV", JOptionPane.ERROR_MESSAGE );
      return;
    }
    JOptionPane.showMessageDialog( parent, "Salvo em " + f.getAbsolutePath(), "CSV", JOptionPane.INFORMATION_MESSAGE );
  }
  private static String csv( double sec ) {
    return Double.isNaN( sec ) ? "" : String.format( Locale.US, "%.4f", 1000.0 * sec );
  }

  //--------------------------------------------------------------------
  // Header detection
  //--------------------------------------------------------------------
  /** Candidate names of the node header, in order of preference (case-insensitive) */
  static final String[] NODE_HDRS = { "rcv", "station", "rcv_station", "rec_station", "node", "node_id", "receiver" };
  /** Candidate names of the sensor/component header */
  static final String[] SENSOR_HDRS = { "sensor", "sensor_id", "comp", "component", "chan", "channel" };

  /** Index of the first header matching `preferred` (last used) or one of `candidates`, exact
   *  name first (case-insensitive, trimmed), then any header NAME containing the candidate; -1 if none */
  static int detectHeader( csHeaderDef[] defs, String preferred, String[] candidates ) {
    List<String> all = new ArrayList<>();
    if( preferred != null ) all.add( preferred );
    for( String c : candidates ) all.add( c );
    for( String c : all ) {
      for( int i = 0; i < defs.length; i++ ) {
        if( defs[i].name != null && defs[i].name.trim().equalsIgnoreCase( c ) ) return i;
      }
    }
    for( String c : candidates ) {
      for( int i = 0; i < defs.length; i++ ) {
        if( defs[i].name != null && defs[i].name.toLowerCase( Locale.US ).contains( c ) ) return i;
      }
    }
    return -1;
  }
  /** Distinct sensor values in trace order (with no sensor header: position of the trace inside its node, 1..m) */
  static List<String> sensorValues( csISeismicTraceBuffer buf, int iNode, int iSensor ) {
    Set<String> vals = new LinkedHashSet<>();
    Map<Integer,Integer> count = new LinkedHashMap<>();
    for( int i = 0; i < buf.numTraces(); i++ ) {
      csHeader[] h = buf.headerValues( i );
      if( iSensor >= 0 ) vals.add( hdrText( h[iSensor] ) );
      else vals.add( String.valueOf( count.merge( h[iNode].intValue(), 1, Integer::sum ) ) );
    }
    return new ArrayList<>( vals );
  }
  static int countNodes( csISeismicTraceBuffer buf, int iNode ) {
    Set<Integer> s = new LinkedHashSet<>();
    for( int i = 0; i < buf.numTraces(); i++ ) s.add( buf.headerValues( i )[iNode].intValue() );
    return s.size();
  }
  /** Header value as text: strings as they are, numbers without ".0" when integral */
  static String hdrText( csHeader h ) {
    if( h == null ) return "";
    Object v = h.value();
    if( v == null ) return "";
    if( v instanceof Number ) {
      double d = ( (Number)v ).doubleValue();
      if( d == Math.rint( d ) && Math.abs( d ) < 1.0e15 ) return String.valueOf( (long)d );
      return String.valueOf( d );
    }
    return v.toString().trim();
  }
  /** "name = value" of the first `max` headers of one trace */
  static String headerList( csHeaderDef[] defs, csHeader[] h, int max ) {
    StringBuilder sb = new StringBuilder();
    for( int i = 0; i < defs.length && i < max; i++ ) {
      sb.append( String.format( Locale.US, "  %-18s = %s%n", defs[i].name, ( h != null && i < h.length ) ? hdrText( h[i] ) : "" ) );
    }
    if( defs.length > max ) sb.append( String.format( Locale.US, "  ... (+%d)%n", defs.length - max ) );
    return sb.toString();
  }
  /** "Headers..." dialog: every header of the 1st trace, with node/sensor choosers. @return true if changed */
  private static boolean chooseHeaders( csIPluginContext ctx, csHeaderDef[] defs, csISeismicTraceBuffer buf, int[] sel ) {
    String[] names = new String[defs.length];
    for( int i = 0; i < defs.length; i++ ) names[i] = defs[i].name;
    String[] sNames = new String[names.length + 1];
    sNames[0] = "(nenhum: ordem do traço no node)";
    System.arraycopy( names, 0, sNames, 1, names.length );
    JComboBox<String> cNode = new JComboBox<>( names );
    cNode.setSelectedIndex( sel[0] );
    JComboBox<String> cSens = new JComboBox<>( sNames );
    cSens.setSelectedIndex( sel[1] + 1 );
    cNode.setMaximumRowCount( 25 );
    cSens.setMaximumRowCount( 25 );
    JTextArea list = new JTextArea( headerList( defs, buf.headerValues( 0 ), defs.length ), 18, 46 );
    list.setEditable( false );
    list.setFont( new Font( Font.MONOSPACED, Font.PLAIN, 12 ) );
    list.setCaretPosition( 0 );
    JPanel choose = new JPanel( new GridLayout( 0, 2, 6, 4 ) );
    choose.add( new JLabel( "Header do node:" ) );   choose.add( cNode );
    choose.add( new JLabel( "Header do sensor:" ) ); choose.add( cSens );
    JPanel p = new JPanel( new BorderLayout( 4, 8 ) );
    p.add( new JLabel( String.format( Locale.US, "%d headers no painel (valores do 1º traço):", defs.length ) ), BorderLayout.NORTH );
    p.add( new JScrollPane( list ), BorderLayout.CENTER );
    p.add( choose, BorderLayout.SOUTH );
    if( JOptionPane.showConfirmDialog( ctx.getParentFrame(), p, "Headers", JOptionPane.OK_CANCEL_OPTION,
        JOptionPane.PLAIN_MESSAGE ) != JOptionPane.OK_OPTION ) return false;
    sel[0] = cNode.getSelectedIndex();
    sel[1] = cSens.getSelectedIndex() - 1;
    return true;
  }

  //--------------------------------------------------------------------
  private static boolean sameValue( String a, String b ) {
    if( a.equals( b ) ) return true;
    try { return Double.parseDouble( a ) == Double.parseDouble( b ); }
    catch( NumberFormatException e ) { return false; }
  }
  private static JTextField num( double v ) {
    return new JTextField( ( v == Math.rint( v ) && Math.abs( v ) < 1e9 ) ? String.valueOf( (long)v ) : String.valueOf( v ) );
  }
  private static double dbl( JTextField f ) {
    return Double.parseDouble( f.getText().trim().replace( ',', '.' ) );
  }
  private static csClockXcorr.Params copy( csClockXcorr.Params s ) {
    csClockXcorr.Params p = new csClockXcorr.Params();
    p.fmin = s.fmin; p.fmax = s.fmax; p.winSec = s.winSec; p.maxLagSec = s.maxLagSec; p.minLagSec = s.minLagSec;
    p.periodMin = s.periodMin; p.oneBit = s.oneBit; p.whiten = s.whiten; p.neighbors = s.neighbors;
    p.minQuality = s.minQuality; p.refNode = s.refNode;
    return p;
  }
}
