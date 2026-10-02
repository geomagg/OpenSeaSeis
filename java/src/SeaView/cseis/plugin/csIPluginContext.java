/* SeaView plugin context: M. A. Gallotti Guimaraes (COPPE-UFRJ), 2026 */

package cseis.plugin;

import cseis.processing.csIProcessing;
import cseis.seis.csHeaderDef;
import cseis.seis.csISeismicReader;
import cseis.seis.csISeismicTraceBuffer;
import javax.swing.JFrame;

/**
 * What a plugin can see and do in SeaView: read the traces of the active pane,
 * open a new pane with computed data, or attach an on-screen processing step.
 */
public interface csIPluginContext {
  /** @return SeaView main window, to be used as parent of dialogs */
  public JFrame getParentFrame();
  /** @return Traces currently loaded in the active pane (original, unprocessed samples) */
  public csISeismicTraceBuffer getTraceBuffer();
  /** @return Trace header definitions of the active pane */
  public csHeaderDef[] getHeaderDef();
  /** @return Index of trace header with the given name, or -1 if it does not exist */
  public int getHeaderIndex( String name );
  /** @return Sample interval of the active pane [ms], [m] or [Hz] */
  public float getSampleInt();
  /** @return Number of samples per trace */
  public int getNumSamples();
  /** @return Vertical domain: cseis.general.csUnits.DOMAIN_TIME/DOMAIN_DEPTH/DOMAIN_FREQ */
  public int getVerticalDomain();
  /** @return Title (file name) of the active pane */
  public String getTitle();
  /**
   * Open computed data in a new SeaView pane
   * @param reader Data, e.g. a cseis.jni.csVirtualSeismicReader filled by the plugin
   * @param title  Title of the new pane
   */
  public void openNewPane( csISeismicReader reader, String title );
  /**
   * Attach an on-screen processing step to the active pane and open its parameter dialog
   * (as the AGC and filter steps of menu 'Active pane / Processing').
   */
  public void addProcessingStep( csIProcessing proc );
}
