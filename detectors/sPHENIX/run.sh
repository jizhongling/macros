#!/bin/bash

PROC=$1
NEV=$2

mkdir -p $SPIN/data/sphenix/output $SPIN/data/sphenix/histos
root -l -b -q Fun4All_G4_sPHENIX.C\($PROC,$NEV\)
root -l -b -q ReadHitFile.C\($PROC,$NEV\)
