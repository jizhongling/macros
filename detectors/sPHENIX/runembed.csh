#!/bin/tcsh

# Input parameters
set proc = $1
set nevents = $2
set runno = 62
@ index = $proc * $nevents / 360
@ skip = $proc * $nevents % 360

# Construct the G4Hits DST files to access. These are MinBias 50 kHz pile up AuAu events
# Lustre location: /sphenix/lustre01/sphnxpro/mdc2/shijing_hepmc/fm_0_20/trkrhit
set strembed0=`printf "DST_TRUTH_sHijing_0_20fm_50kHz_bkg_0_20fm-%010d-%05d.root" $runno $index`
# Lustre location: /sphenix/lustre01/sphnxpro/mdc2/shijing_hepmc/fm_0_20/pileup
set strembed1=`printf "DST_TRKR_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-%010d-%05d.root" $runno $index`

# Run my Fun4AllMacro which locates the G4Hits files by FROG
mkdir -p $SPIN/data/sphenix/output $SPIN/data/sphenix/histos
root -l -b -q 'Fun4All_G4_sPHENIX.C('$proc', '$nevents', "'$strembed0'", "'$strembed1'", '$skip')'
