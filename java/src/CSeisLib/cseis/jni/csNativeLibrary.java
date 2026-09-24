package cseis.jni;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;

/**
 * Carrega a biblioteca nativa csJNIlib (leitores SEG-Y/SEG-D/RSF/CSeis, FFT, filtros).
 * <p>
 * Ordem de busca:
 * <ol>
 * <li>Propriedade {@code cseis.native.lib}: caminho completo do arquivo da biblioteca.</li>
 * <li>{@code java.library.path} (todas as entradas, via {@link System#loadLibrary}).</li>
 * <li>Cópia embutida no jar em {@code /native/<os>-<arch>/}, extraída para um diretório temporário.</li>
 * </ol>
 * Chamadas repetidas são seguras: a biblioteca é carregada uma única vez.
 */
public final class csNativeLibrary {
  public static final String LIB_NAME = "csJNIlib";

  private static boolean loaded = false;
  private static String loadedFrom = null;

  private csNativeLibrary() {}

  /** Carrega a biblioteca; lança UnsatisfiedLinkError com mensagem detalhada se falhar. */
  public static synchronized void load() {
    if( loaded ) return;
    StringBuilder errors = new StringBuilder();

    String explicit = System.getProperty( "cseis.native.lib" );
    if( explicit != null && !explicit.isEmpty() ) {
      try {
        System.load( new File( explicit ).getAbsolutePath() );
        markLoaded( explicit );
        return;
      }
      catch( UnsatisfiedLinkError e ) {
        errors.append( "cseis.native.lib=" ).append( explicit ).append( ": " ).append( e.getMessage() ).append( '\n' );
      }
    }

    try {
      System.loadLibrary( LIB_NAME );
      markLoaded( "java.library.path" );
      return;
    }
    catch( UnsatisfiedLinkError e ) {
      errors.append( "java.library.path: " ).append( e.getMessage() ).append( '\n' );
    }

    String resource = "/native/" + platformId() + "/" + System.mapLibraryName( LIB_NAME );
    try( InputStream in = csNativeLibrary.class.getResourceAsStream( resource ) ) {
      if( in == null ) {
        errors.append( "jar: recurso " ).append( resource ).append( " não encontrado (plataforma sem binário embutido)\n" );
      }
      else {
        Path dir = Files.createTempDirectory( "seaview-native" );
        Path lib = dir.resolve( System.mapLibraryName( LIB_NAME ) );
        Files.copy( in, lib, StandardCopyOption.REPLACE_EXISTING );
        lib.toFile().deleteOnExit();
        dir.toFile().deleteOnExit();
        System.load( lib.toAbsolutePath().toString() );
        markLoaded( "jar (" + resource + ")" );
        return;
      }
    }
    catch( IOException | UnsatisfiedLinkError e ) {
      errors.append( "jar: " ).append( e.getMessage() ).append( '\n' );
    }

    throw new UnsatisfiedLinkError( "Não foi possível carregar a biblioteca nativa '" + LIB_NAME + "'.\n" + errors );
  }

  /** Tenta carregar; retorna false em vez de lançar exceção. */
  public static boolean tryLoad() {
    try {
      load();
      return true;
    }
    catch( UnsatisfiedLinkError e ) {
      System.err.println( e.getMessage() );
      return false;
    }
  }

  public static synchronized boolean isLoaded() { return loaded; }
  public static synchronized String loadedFrom() { return loadedFrom; }

  /** Identificador da plataforma, ex.: linux-x86_64, windows-x86_64, macos-aarch64. */
  public static String platformId() {
    String os = System.getProperty( "os.name", "" ).toLowerCase();
    String arch = System.getProperty( "os.arch", "" ).toLowerCase();
    String o = os.startsWith( "windows" ) ? "windows" : os.startsWith( "mac" ) ? "macos" : "linux";
    String a = ( arch.equals( "amd64" ) || arch.equals( "x86_64" ) ) ? "x86_64"
             : ( arch.equals( "aarch64" ) || arch.equals( "arm64" ) ) ? "aarch64" : arch;
    return o + "-" + a;
  }

  private static void markLoaded( String from ) {
    loaded = true;
    loadedFrom = from;
  }
}
