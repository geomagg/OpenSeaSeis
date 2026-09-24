package cseis.seaview.plugin;

/**
 * Ponto de extensão do SeaView.
 * <p>
 * Para criar uma nova funcionalidade:
 * <ol>
 * <li>Crie a classe em {@code java/src/SeaView/cseis/plugins/}.</li>
 * <li>Adicione o nome completo da classe em {@code java/src/SeaView/cseis/resources/seaview_plugins.txt}.</li>
 * <li>Rode {@code ./make_java.sh} e copie os jars para a pasta lib da instalação.</li>
 * </ol>
 * O método {@link #install} roda na thread do Swing, uma vez, ao abrir o SeaView.
 */
public interface csSeaViewPlugin {

  /** Nome curto, usado no menu e nas mensagens de erro. */
  String getName();

  /** Registra menus/ações. Use {@link csPluginContext#addMenuItem} para aparecer no menu "Plugins". */
  void install( csPluginContext context );
}
