#!/bin/bash

cd ~/freesurfer_repo_root/freesurfer/mris_fix_topology_test
rm -rf results* subjects*
tar xzf 000_results_subjects.tgz

export FREESURFER_HOME=/home/rheticus/freesurfer_repo_root/freesurfer/distribution
export FS_LICENSE=~/license.txt
export SUBJECTS_DIR=subjects
export FREESURFER_REPLACEMENT_FOR_CREATION_TIME_STRING="Day Mon 00 00:00:00 0000"

for numThreads in 4 2 1 ; do

    rm -rf   /tmp/ROMP_statsFiles/
    mkdir -p /tmp/ROMP_statsFiles/

    export OMP_NUM_THREADS=$numThreads

    ../mris_fix_topology/mris_fix_topology \
            -niters 1 \
            -mgz \
            -sphere qsphere.nofix \
            -ga -seed 1234 bert lh > results/output.txt

    strings ./subjects/bert/surf/lh.orig > results/lh.orig.strings

    mv results{,$numThreads}
    mv subjects{,$numThreads}

done
