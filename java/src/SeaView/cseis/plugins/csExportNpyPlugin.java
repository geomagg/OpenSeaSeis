package cseis.plugins;

import cseis.seaview.csSeisPaneBundle;
import cseis.seaview.plugin.csPluginContext;
import cseis.seaview.plugin.csSeaViewPlugin;
import cseis.seis.csHeader;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicTraceBuffer;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.util.Locale;

/**
 * Exemplo de plugin: exporta os traços do painel ativo para NumPy.
 * <ul>
 * <li>{@code nome.npy}: float32, shape (nTraços, nAmostras) — {@code np.load("nome.npy")}</li>
 * <li>{@code nome_headers.csv}: um traço por linha, uma coluna por cabeçalho — {@code pandas.read_csv}</li>
 * </ul>
 */
public class csExportNpyPlugin implements csSeaViewPlugin {
  private File lastDir = null;

  @Override
  public String getName() { return "Exportar para NumPy"; }

  @Override
  public void install( csPluginContext ctx ) {
    ctx.addBundleMenuItem( "Exportar traços do painel (.npy + cabeçalhos .csv)...", "ctrl shift E", bundle -> export( ctx, bundle ) );
  }

  private void export( csPluginContext ctx, csSeisPaneBundle bundle ) {
    JFileChooser fc = new JFileChooser( lastDir != null ? lastDir : new File( bundle.getPath() ) );
    fc.setFileFilter( new FileNameExtensionFilter( "NumPy array (*.npy)", "npy" ) );
    String base = new File( bundle.getFilenamePath() ).getName().replaceFirst( "\\.[^.]+$", "" );
    fc.setSelectedFile( new File( fc.getCurrentDirectory(), base + ".npy" ) );
    if( fc.showSaveDialog( ctx.getSeaView() ) != JFileChooser.APPROVE_OPTION ) return;

    File npy = fc.getSelectedFile();
    if( !npy.getName().toLowerCase().endsWith( ".npy" ) ) npy = new File( npy.getPath() + ".npy" );
    lastDir = npy.getParentFile();
    File csv = new File( npy.getParentFile(), npy.getName().replaceFirst( "\\.npy$", "" ) + "_headers.csv" );

    csISeismicTraceBuffer buf = bundle.getTraceBuffer();
    try {
      writeNpy( npy, buf );
      writeHeadersCsv( csv, buf, bundle.getTraceHeaderDefs() );
    }
    catch( IOException e ) {
      ctx.error( "Erro ao exportar:\n" + e.getMessage() );
      return;
    }
    ctx.setStatus( "Exportado: " + npy.getName() + " (" + buf.numTraces() + " x " + buf.numSamples() + ") e " + csv.getName() );
    ctx.info( String.format( Locale.US,
        "Exportados %d traços x %d amostras (dt = %g ms).%n%n%s%n%s%n%nPython:%n  import numpy as np, pandas as pd%n  d = np.load(r\"%s\")%n  h = pd.read_csv(r\"%s\")",
        buf.numTraces(), buf.numSamples(), bundle.getSampleInt(), npy.getPath(), csv.getPath(), npy.getPath(), csv.getPath() ) );
  }

  /** Escreve array float32 little-endian no formato .npy versão 1.0. */
  static void writeNpy( File file, csISeismicTraceBuffer buf ) throws IOException {
    int ntr = buf.numTraces();
    int ns = buf.numSamples();
    String dict = "{'descr': '<f4', 'fortran_order': False, 'shape': (" + ntr + ", " + ns + "), }";
    int preamble = 6 + 2 + 2;                              // magic + versão + tamanho do header
    int total = preamble + dict.length() + 1;              // +1 do '\n'
    int pad = ( 64 - total % 64 ) % 64;
    StringBuilder header = new StringBuilder( dict );
    for( int i = 0; i < pad; i++ ) header.append( ' ' );
    header.append( '\n' );

    try( OutputStream os = new BufferedOutputStream( new FileOutputStream( file ), 1 << 20 ) ) {
      os.write( new byte[]{ (byte)0x93, 'N', 'U', 'M', 'P', 'Y', 1, 0 } );
      int hlen = header.length();
      os.write( hlen & 0xff );
      os.write( ( hlen >> 8 ) & 0xff );
      os.write( header.toString().getBytes( StandardCharsets.US_ASCII ) );
      ByteBuffer bb = ByteBuffer.allocate( ns * 4 ).order( ByteOrder.LITTLE_ENDIAN );
      for( int i = 0; i < ntr; i++ ) {
        bb.clear();
        float[] s = buf.samples( i );
        for( int k = 0; k < ns; k++ ) bb.putFloat( k < s.length ? s[k] : 0f );
        os.write( bb.array() );
      }
    }
  }

  static void writeHeadersCsv( File file, csISeismicTraceBuffer buf, csHeaderDef[] defs ) throws IOException {
    try( PrintWriter pw = new PrintWriter( new java.io.OutputStreamWriter( new FileOutputStream( file ), StandardCharsets.UTF_8 ) ) ) {
      StringBuilder line = new StringBuilder( "trace_number" );
      for( csHeaderDef d : defs ) line.append( ',' ).append( d.name );
      pw.println( line );
      for( int i = 0; i < buf.numTraces(); i++ ) {
        line.setLength( 0 );
        line.append( buf.originalTraceNumber( i ) );
        csHeader[] h = buf.headerValues( i );
        for( int k = 0; k < defs.length; k++ ) {
          line.append( ',' );
          if( h != null && k < h.length && h[k] != null ) {
            Object v = h[k].value();
            String s = String.valueOf( v );
            if( v instanceof String && ( s.contains( "," ) || s.contains( "\"" ) ) ) s = '"' + s.replace( "\"", "\"\"" ) + '"';
            line.append( s );
          }
        }
        pw.println( line );
      }
    }
  }
}
