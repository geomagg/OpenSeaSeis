#!/bin/bash
# modela o no n (1..17): fonte no no (coluna 129+16*(n-1), 990 m), receptores a 10 m
cd "$(dirname "$0")"
n=$1; x=$((129 + (n-1)*16))
seaseis -D X=$x -D N=$n -f fd.flow -o fd$n.log > /dev/null 2>&1
echo "no $n (coluna $x): rc=$?"
