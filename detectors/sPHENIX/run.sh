#!/bin/bash

PROC=$1
inFile=$SPIN/data/sphenix/output/G4sPHENIX_g4svtx_eval-$PROC.root
if [ ! -f $inFile ] ; then
  exit 0
fi

mkdir -p $SPIN/data/sphenix/output $SPIN/data/sphenix/histos
for i in {0..4} ; do
  ./AnaHitFile $inFile $SPIN/data/sphenix/histos/training-$PROC $i &
done
wait
