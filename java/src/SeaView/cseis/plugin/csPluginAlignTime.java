/* SeaView plugin: absolute-time alignment. M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import javax.swing.JOptionPane;

/**
 * Plugin: align all traces of the active pane in absolute time (common start = latest start,
 * optional cut at the earliest end), as an on-screen processing step -- see csProcessingAlignTime.
 */
public class csPluginAlignTime implements csISeaViewPlugin {
  @Override
  public String getName() {
    return "Alinhar traços no tempo absoluto...";
  }
  @Override
  public String getDescription() {
    return "Todos os traços passam a começar no mesmo instante (o início mais tardio) e terminar no fim mais cedo; usa time_samp1 ou time_day/hour/min/sec";
  }
  @Override
  public void run( csIPluginContext ctx ) {
    if( ctx.getTraceBuffer() == null || ctx.getTraceBuffer().numTraces() == 0 ) {
      JOptionPane.showMessageDialog( ctx.getParentFrame(), "Sem dados no painel ativo", getName(), JOptionPane.WARNING_MESSAGE );
      return;
    }
    ctx.addProcessingStep( new csProcessingAlignTime( ctx.getHeaderDef(), ctx.getSampleInt(), ctx.getTraceBuffer() ) );
  }
}
