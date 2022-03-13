#!/bin/tcsh

# input parameters
set runno = $1
set nevents = $2
@ fileno = $runno * $nevents / 500
@ skip = $runno * $nevents % 500
if ($fileno >= 10) then
  @ fileno += 3300 - 10
endif

# Construct the G4Hits files to access. These are AuAu events
set strembed=`printf "G4Hits_sHijing_0_20fm-0000000003-%05d.root" $fileno`

# Construct the G4Hits DST files to access. These are MinBias 50 kHz pile up AuAu events
#set strembed0=`printf "DST_TRUTH_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000003-%05d.root" $fileno`
#set strembed1=`printf "DST_TRKR_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000003-%05d.root" $fileno`

# Run my Fun4AllMacro which locates the G4Hits files by FROG
mkdir -p $SPIN/data/sphenix/output $SPIN/data/sphenix/histos
root -l -b -q 'Fun4All_G4_sPHENIX.C('$runno', '$nevents', "'$strembed'", "", '$skip')'
