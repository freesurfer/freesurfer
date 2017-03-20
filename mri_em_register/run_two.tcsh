#!/bin/tcsh -f

cd testdata

echo "Running two mri_em_register"

$PERF ../mri_em_register \
        -uns 3 \
        -mask brainmask.mgz \
        nu.mgz \
        ../../distribution/average/RB_all_2016-05-10.vc700.gca \
        talairach.run_two_0.lta >& ../run_two_${OMP_NUM_THREADS}_0.log &

$PERF ../mri_em_register \
        -uns 3 \
        -mask brainmask.mgz \
        nu.mgz \
        ../../distribution/average/RB_all_2016-05-10.vc700.gca \
        talairach.run_two_1.lta >& ../run_two_${OMP_NUM_THREADS}_1.log &

wait
