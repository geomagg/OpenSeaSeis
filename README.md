## Fork notes (geomagg)

- SeaView **Plugins** menu: add a class in `java/src/SeaView/cseis/plugins/`, list it in
  `java/src/SeaView/cseis/resources/seaview_plugins.txt`, then run `./make_java.sh` and `./make_seaseis.sh`.
  Examples: panel statistics (Ctrl+Shift+I) and NumPy export (Ctrl+Shift+E).
- `csJNIlib` is loaded through `cseis.jni.csNativeLibrary`, which searches every `java.library.path` entry.
- `src/cs/jni`: fixed `GetMethodID` failure checks, which modern GCC rejects.
- `make_java.sh`: compiles with `-encoding UTF-8`.

# OpenSeaSeis

OpenSeaSeis is based on the package SeaSeis that was created by Bjorn Olofsen in 2006. The package is now sole property
of the Colorado School of Mines.


See subdirectory 'doc' for more information

README_OPENSEASEIS_V3.00 - General information about this release, and short ins
tall notes
README_SEISMIC_UNIX  - Seaseis and Seismic Unix

Seaview_tutorial_v1.61.pdf - Tutorial from 2011 (Seaview is a seismic 2D viewer)
 ...some of Seaview's features have changed since then but most are still the sa
me.

EAGE 2012 Open-software Workshop abstract and poster
  eage2012_open_source_workshop_Seaseis.pdf
  eage2012_open_source_workshop_Seaseis_poster.pdf

