#!/bin/tcsh -f

cd testdata

echo "Running single mri_em_register"

$PERF ../mri_em_register \
        -uns 3 \
        -mask brainmask.mgz \
        nu.mgz \
        ../../distribution/average/RB_all_2016-05-10.vc700.gca \
        talairach.run_one.lta >& ../run_one_${OMP_NUM_THREADS}.log &

wait

