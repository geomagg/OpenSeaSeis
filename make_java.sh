#!/bin/bash

# Usage:   make_java.sh [option]
# Options: clean

JAVAC=javac

source src/make/linux/set_environment.sh

export CLASS_CSEISLIB=${CSEISDIR}/class_CSeisLib
export CLASS_SEAVIEW=${CSEISDIR}/class_SeaView
export CLASS_XCSEIS=${CSEISDIR}/class_XCSeis

if [ $# -eq 1 ]; then
  if [ "$1" == "clean" ]; then
    echo "Cleaning directory ${CLASS_CSEISLIB}..."
    rm -r ${CLASS_CSEISLIB}
    echo "Cleaning directory ${CLASS_SEAVIEW}..."
    rm -r ${CLASS_SEAVIEW}
    echo "Cleaning directory ${CLASS_XCSEIS}..."
    rm -r ${CLASS_XCSEIS}
  else
    echo "Option not recognized: $1"
  fi
  exit -1
fi

mkdir -p ${CLASS_CSEISLIB}
mkdir -p ${CLASS_SEAVIEW}
mkdir -p ${CLASS_XCSEIS}

echo ${CLASS_CSEISLIB}

echo "Compile CSeisLib..."
$JAVAC -encoding UTF-8 java/src/CSeisLib/cseis/*/*.java -sourcepath java/src/CSeisLib/ -d ${CLASS_CSEISLIB}
echo "Compile SeaView..."
$JAVAC -encoding UTF-8 java/src/SeaView/cseis/*/*.java -sourcepath java/src/SeaView/ -d ${CLASS_SEAVIEW} -cp ${CLASS_CSEISLIB}
echo "Compile XCSeis..."
$JAVAC -encoding UTF-8 java/src/XCSeis/cseis/*/*.java -sourcepath java/src/XCSeis/ -d ${CLASS_XCSEIS} -cp ${CLASS_CSEISLIB}:${CLASS_SEAVIEW}

mkdir -p ${CLASS_CSEISLIB}/cseis/resources
mkdir -p ${CLASS_SEAVIEW}/cseis/resources
mkdir -p ${CLASS_XCSEIS}/cseis/resources

cp java/src/CSeisLib/cseis/resources/* ${CLASS_CSEISLIB}/cseis/resources/
cp java/src/SeaView/cseis/resources/* ${CLASS_SEAVIEW}/cseis/resources/
#cp java/src/XCSeis/cseis/resources/* ${CLASS_XCSEIS}/cseis/resources/

echo "Building jar files..."

jar cf java/jar/CSeisLib.jar -C ${CLASS_CSEISLIB} cseis
jar cf java/jar/SeaView.jar -C ${CLASS_SEAVIEW} cseis
jar cf java/jar/XCSeis.jar -C ${CLASS_XCSEIS} cseis

echo "Done."
