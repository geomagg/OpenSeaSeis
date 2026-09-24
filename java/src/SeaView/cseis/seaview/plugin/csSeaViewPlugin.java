package cseis.seaview.plugin;

/**
 * Ponto de extensão do SeaView.
 * <p>
 * Para criar uma nova funcionalidade:
 * <ol>
 * <li>Implemente esta interface (ex.: no módulo {@code plugins/}).</li>
 * <li>Registre a classe em {@code META-INF/services/cseis.seaview.plugin.csSeaViewPlugin}
 *     (uma classe por linha).</li>
 * <li>Recompile. Plugins do módulo {@code plugins/} já vão dentro do executável;
 *     jars externos podem ser colocados na pasta {@code plugins/} ao lado do executável
 *     ou em {@code ~/.seaview/plugins}.</li>
 * </ol>
 * O método {@link #install} roda na thread do Swing, uma vez, ao abrir o SeaView.
 */
public interface csSeaViewPlugin {

  /** Nome curto, usado no menu e nas mensagens de erro. */
  String getName();

  /** Registra menus/ações. Use {@link csPluginContext#addMenuItem} para aparecer no menu "Plugins". */
  void install( csPluginContext context );
}
