#!/bin/tcsh

# Input parameters
set runno = $1
set nevents = $2
@ fileno = $runno * $nevents / 400
@ skip = $runno * $nevents % 400
if ($fileno >= 0) then
  @ fileno += 3300
endif

# Construct the G4Hits DST files to access. These are MinBias 50 kHz pile up AuAu events
set strembed0=`printf "DST_TRUTH_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000003-%05d.root" $fileno`
set strembed1=`printf "DST_TRKR_HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000003-%05d.root" $fileno`

# Run my Fun4AllMacro which locates the G4Hits files by FROG
mkdir -p $SPIN/data/sphenix/output $SPIN/data/sphenix/histos
root -l -b -q 'Fun4All_G4_sPHENIX.C('$runno', '$nevents', "'$strembed0'", "'$strembed1'", '$skip')'
