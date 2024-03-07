#!/bin/bash

proc=$1
tree_dir=$TGHF/data/sphenix/output
hist_dir=$SPIN/data/sphenix/histos
#mkdir -p $tree_dir $hist_dir

for i in {0..1} ; do
  #./AnaToyClusters 10000 $hist_dir/training-$proc $i &
  #./AnaHitFile $tree_dir/G4sPHENIX_g4svtx_eval-$proc.root $hist_dir/training-$proc $i 5 &
  #./AnaQA $tree_dir/G4sPHENIX_g4svtx_eval- $(( proc*20+i*10 )) $(( proc*20+i*10+10 )) $hist_dir/qa-$proc-$i.root &
  #./AnaPull $tree_dir/G4sPHENIX_g4svtx_eval- $(( proc*20+i*10 )) $(( proc*20+i*10+10 )) $histo_dir/pull-$proc-$i.root &
  #./AnaHFElectron $tree_dir $(( proc*20+i*10 )) $(( proc*20+i*10+10 )) $hist_dir/hf-electron $proc-$i &
  ./AnaEEMass $tree_dir $(( proc*20+i*10 )) $(( proc*20+i*10+10 )) $hist_dir/ee-minv $proc-$i &
done
wait
