#!/bin/bash
# Build the example SeaView plugin and install it in the plugin directory.
# Usage: ./make_plugin.sh            (needs CSEISDIR, e.g. from ~/.bashrc)
# The jar goes to $CSEISDIR/lib/plugins; SeaView shows the plugin in menu 'Plugins' at the next start.
set -e
cd "$(dirname "$0")"
LIB=${CSEISDIR:?set CSEISDIR}/lib
NAME=$(basename "$(pwd)")
rm -rf classes && mkdir classes
javac -cp "$LIB/CSeisLib.jar:$LIB/SeaView.jar" -d classes $(find src -name "*.java")
cp -r META-INF classes/
mkdir -p "$LIB/plugins"
jar cf "$LIB/plugins/seaview_plugin_$NAME.jar" -C classes .
rm -rf classes
echo "Installed $LIB/plugins/seaview_plugin_$NAME.jar"
