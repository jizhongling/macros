#!/bin/tcsh

# Input parameters
set runno = $1
set nevents = $2
set datadir = $_CONDOR_SCRATCH_DIR/proc$runno

# Construct the G4Hits DST files to access. These are MinBias 50 kHz pile up AuAu events
set strembed0=`printf "DST_TRUTH_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000002-%05d.root" $runno`
set strembed1=`printf "DST_TRKR_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000002-%05d.root" $runno`

# Copy the sPHENIX get input files perl script to the condor scratch directory to download the G4Hits files from dCache
mkdir -p $datadir
cp getinputfiles.pl $datadir

pushd $datadir
chmod +x getinputfiles.pl

# Put the DSTs into a filelist.txt
echo $strembed0 > filelist.txt
echo $strembed1 >> filelist.txt

# Run the perl script and download the files from dcache
getinputfiles.pl --dcache --filelist filelist.txt
ls -lhrt
popd

# Run my Fun4AllMacro pointing to the downloaded G4Hits files in the condor scratch directory
mkdir -p $SPIN/data/sphenix/output $SPIN/data/sphenix/histos
root -l -b -q 'Fun4All_G4_sPHENIX.C('$runno', '$nevents', "'$datadir/$strembed0'", "'$datadir/$strembed1'")'
