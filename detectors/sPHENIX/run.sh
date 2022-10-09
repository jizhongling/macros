#!/bin/bash

PROC=$1
PREFIX=$SPIN/data/sphenix/output/G4sPHENIX_g4svtx_eval-

mkdir -p $SPIN/data/sphenix/output $SPIN/data/sphenix/histos
for i in {0..1} ; do
  #./AnaToyClusters 10000 $SPIN/data/sphenix/histos/training-$PROC $i &
  ./AnaHitFile $PREFIX$PROC.root $SPIN/data/sphenix/histos/training-$PROC $i &
  #./AnaQA $PREFIX $(( PROC*20+i*10 )) $(( PROC*20+i*10+10 )) $SPIN/data/sphenix/histos/qa-$PROC-$i.root &
  #./AnaPull $PREFIX $(( PROC*20+i*10 )) $(( PROC*20+i*10+10 )) $SPIN/data/sphenix/histos/pull-$PROC-$i.root &
done
wait
