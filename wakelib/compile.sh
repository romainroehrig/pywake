#!/bin/sh

set -ex

solib='wake'
routines='parkind1 yomcst sucst yomwake suwake yomphy2 suphy2 wake'

for rr in $routines
do
  gfortran -fbacktrace -fPIC -c $rr.F90 -o $rr.o
done

gfortran -fbacktrace -fPIC -shared -O2 *.o -o $solib.so

rm -f *.o *.mod
