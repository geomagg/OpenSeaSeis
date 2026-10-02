SeaView plugins
===============

SeaView has a menu 'Plugins'. Each plugin acts on the active data pane.

Built-in plugins (source in java/src/SeaView/cseis/plugin):
  F-X spectrum...          amplitude spectrum of each trace in a time window, in a new pane
  HMO / NMO (velocity)...  hyperbolic moveout with a chosen velocity, applied on screen

External plugins are jar files placed in one of these directories:
  $CSEISDIR/lib/plugins        (next to SeaView.jar)
  ~/.seaview/plugins
  the directory in the environment variable SEAVIEW_PLUGINS
Menu 'Plugins / About plugins...' lists the directories and the plugins found.

Writing a plugin
----------------
1. Copy the directory 'example' to a new name.
2. Write a class implementing cseis.plugin.csISeaViewPlugin (getName, getDescription, run).
   run() gets a cseis.plugin.csIPluginContext with the traces, headers and sample interval of the
   active pane, and methods to open a new pane with computed data (openNewPane) or to attach an
   on-screen processing step (addProcessingStep, see csProcessingHMO).
3. List the full class name(s) in META-INF/services/cseis.plugin.csISeaViewPlugin.
4. Run ./make_plugin.sh and restart SeaView.
