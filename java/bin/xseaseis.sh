#!/bin/bash
# Launch Job Builder
#

DIR=$(dirname $0)
# Capture the case when script is executed from current directory
LETTER=$(echo $DIR | cut -c 1)
if [ $LETTER == '.' ]; then
   DIR=$(pwd)/$DIR
fi
export LD_LIBRARY_PATH=$DIR/../lib:$LD_LIBRARY_PATH

# Java 17+: allow JNI (csJNIlib) without the "restricted method" warning of Java 24+
NATIVE_ACCESS=""
if java --enable-native-access=ALL-UNNAMED -version >/dev/null 2>&1; then
   NATIVE_ACCESS="--enable-native-access=ALL-UNNAMED"
fi

java ${NATIVE_ACCESS} -cp .:${DIR}/../lib/CSeisLib.jar:${DIR}/../lib/XCSeis.jar:${DIR}/../lib/SeaView.jar -Xmx1524m -Djava.library.path=${DIR}/../lib cseis/xcseis/XCSeis  $*
