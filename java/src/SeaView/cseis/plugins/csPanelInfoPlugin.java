package cseis.plugins;

import cseis.seaview.csSeisPaneBundle;
import cseis.seaview.plugin.csPluginContext;
import cseis.seaview.plugin.csSeaViewPlugin;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicTraceBuffer;

import javax.swing.JDialog;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import java.awt.Font;
import java.util.Locale;

/**
 * Exemplo de plugin: mostra informações e estatísticas de amplitude dos traços carregados no painel ativo.
 */
public class csPanelInfoPlugin implements csSeaViewPlugin {

  @Override
  public String getName() { return "Estatísticas do painel"; }

  @Override
  public void install( csPluginContext ctx ) {
    ctx.addBundleMenuItem( "Estatísticas do painel...", "ctrl shift I", bundle -> show( ctx, bundle ) );
  }

  private void show( csPluginContext ctx, csSeisPaneBundle bundle ) {
    csISeismicTraceBuffer buf = bundle.getTraceBuffer();
    int ntr = buf.numTraces();
    int ns = buf.numSamples();

    double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY, sum = 0, sum2 = 0;
    long n = 0, zeroTraces = 0;
    for( int i = 0; i < ntr; i++ ) {
      float[] s = buf.samples( i );
      boolean allZero = true;
      for( float v : s ) {
        if( v < min ) min = v;
        if( v > max ) max = v;
        sum += v;
        sum2 += (double)v * v;
        if( v != 0f ) allZero = false;
      }
      n += s.length;
      if( allZero ) zeroTraces++;
    }
    double mean = sum / n;
    double rms = Math.sqrt( sum2 / n );

    StringBuilder sb = new StringBuilder();
    sb.append( String.format( Locale.US, "Arquivo:           %s%n", bundle.getFilenamePath() ) );
    sb.append( String.format( Locale.US, "Traços no arquivo: %d%n", bundle.getTotalNumTraces() ) );
    sb.append( String.format( Locale.US, "Traços no painel:  %d  (1º traço original: %d)%n", ntr, buf.originalTraceNumber( 0 ) ) );
    sb.append( String.format( Locale.US, "Amostras/traço:    %d%n", ns ) );
    sb.append( String.format( Locale.US, "Intervalo amostr.: %g %s%n", bundle.getSampleInt(), bundle.isFrequencyDomain() ? "Hz" : "ms" ) );
    sb.append( String.format( Locale.US, "Comprimento:       %g %s%n", bundle.getSampleInt() * ( ns - 1 ), bundle.isFrequencyDomain() ? "Hz" : "ms" ) );
    sb.append( String.format( Locale.US, "%nAmplitude (traços do painel)%n" ) );
    sb.append( String.format( Locale.US, "  mín:   %.6g%n  máx:   %.6g%n  média: %.6g%n  RMS:   %.6g%n", min, max, mean, rms ) );
    sb.append( String.format( Locale.US, "  traços zerados: %d%n", zeroTraces ) );

    csHeaderDef[] hdrs = bundle.getTraceHeaderDefs();
    sb.append( String.format( Locale.US, "%nCabeçalhos de traço (%d):%n", hdrs.length ) );
    for( csHeaderDef h : hdrs ) {
      sb.append( String.format( Locale.US, "  %-20s %s%n", h.name, h.desc == null ? "" : h.desc ) );
    }

    JTextArea text = new JTextArea( sb.toString(), 28, 70 );
    text.setEditable( false );
    text.setFont( new Font( Font.MONOSPACED, Font.PLAIN, 12 ) );
    text.setCaretPosition( 0 );
    JDialog dialog = new JDialog( ctx.getSeaView(), "Estatísticas — " + bundle.getTitle(), false );
    dialog.getContentPane().add( new JScrollPane( text ) );
    dialog.pack();
    dialog.setLocationRelativeTo( ctx.getSeaView() );
    dialog.setVisible( true );
  }
}
