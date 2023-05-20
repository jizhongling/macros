#!/bin/tcsh

# Input parameters
set proc = $1
set nevents = $2
set runno = 6
@ index = $proc * $nevents / 360
@ skip = $proc * $nevents % 360
@ iend = $proc + 1
#set skip = 0

# Output directories
set tree_dir = ${_CONDOR_SCRATCH_DIR}
#set hist_dir = ${_CONDOR_SCRATCH_DIR}
#set tree_dir = $SPIN/data/sphenix/output
set hist_dir = $TGHF/data/sphenix/histos
mkdir -p $tree_dir $hist_dir

# Construct the G4Hits DST files to access. These are MinBias 50 kHz pile up AuAu events
# Lustre location: /sphenix/lustre01/sphnxpro/mdc2/shijing_hepmc/fm_0_20/trkrhit
set strembed0 = `printf "DST_TRUTH_sHijing_0_20fm_50kHz_bkg_0_20fm-%010d-%05d.root" $runno $index`
# Lustre location: /sphenix/lustre01/sphnxpro/mdc2/shijing_hepmc/fm_0_20/pileup
set strembed1 = `printf "DST_TRKR_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-%010d-%05d.root" $runno $index`

# Run the Fun4AllMacro which locates the G4Hits files by FROG
root -l -b -q 'Fun4All_G4_sPHENIX.C('$proc', '$nevents', "'$strembed0'", "'$strembed1'", '$skip', "'$tree_dir'")'

# Run the analysis code
./AnaHitFile $tree_dir/G4sPHENIX_g4svtx_eval-$proc.root $hist_dir/training-$proc 0 $nevents
./AnaQA $tree_dir/G4sPHENIX_g4svtx_eval- $proc $iend $hist_dir/qa-$proc-0.root
