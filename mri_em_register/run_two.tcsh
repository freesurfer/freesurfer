#!/bin/tcsh -f

cd testdata

echo "Running two mri_em_register"

../mri_em_register \
        -uns 3 \
        -mask brainmask.mgz \
        nu.mgz \
        ../../distribution/average/RB_all_2016-05-10.vc700.gca \
        talairach.run_two_0.lta >& ../run_two_0.log &

../mri_em_register \
        -uns 3 \
        -mask brainmask.mgz \
        nu.mgz \
        ../../distribution/average/RB_all_2016-05-10.vc700.gca \
        talairach.run_two_1.lta >& ../run_two_1.log &

wait