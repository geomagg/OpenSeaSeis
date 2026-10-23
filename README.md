## Fork notes (geomagg)

- SeaView **Plugins** menu: add a class in `java/src/SeaView/cseis/plugins/`, list it in
  `java/src/SeaView/cseis/resources/seaview_plugins.txt`, then run `./make_java.sh` and `./make_seaseis.sh`.
  Included: panel info, export to NumPy, and **Spectrogram** (Ctrl+Shift+S, then left-click a trace: STFT with Hann window, dB, linear/log frequency, 8 colormaps, optional detrend and zero-phase Butterworth high/low-pass, PNG export).
  Examples: panel statistics (Ctrl+Shift+I) and NumPy export (Ctrl+Shift+E).
- `csJNIlib` is loaded through `cseis.jni.csNativeLibrary`, which searches every `java.library.path` entry.
- `src/cs/jni`: fixed `GetMethodID` failure checks, which modern GCC rejects.
- `make_java.sh`: compiles with `-encoding UTF-8`.

## SEISPAK modules (SPK_*)

The `SPK_*` modules port programs of **SEISPAK**, the Texaco processing system, into SeaSeis:
binning, gain, deconvolution, filtering, mutes, PZ summation, and 2D and 3D prestack depth
migration of ocean-bottom node data (SPK_XVSPMIG1, SPK_XVSPMIG2, SPK_3DBCPMIG). The
numerical routines are the original Fortran, kept unchanged wherever possible; only the
SEISPAK file I/O was replaced by SeaSeis traces and memory. The changes are marked
`SEASEIS` in the code. Each module was checked against the original program or against
synthetic data.

### In memory of Dwight Sukup

Many of these routines carry the name of **Dwight Sukup** (Texaco E&P Research). The 3D
travel-time routine used by SPK_3DBCPMIG still opens with *"From an original by Dwight Sukup
of Texaco E&P Research"*. In 2014 Dwight rewrote these routines for ocean-bottom node (OBN)
data, together with **Patricia Crawford** and **Marcos Gallotti**, both of the former
Geonunes. The node migrations here (SPK_XVSPMIG1, SPK_XVSPMIG2, SPK_3DBCPMIG) come from
that work. The modules keep the SEISPAK name (the `SPK_` prefix) as a tribute to the Texaco
group that created it. Decades after SEISPAK began, the code can still be followed line by line, and it
still produces correct images. This fork is dedicated to his memory and to everyone who
worked with him on SEISPAK.

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

