#!/usr/bin/env bash
source "$(dirname $0)/../../../test.sh"

for threads in 1 8; do
    export OMP_NUM_THREADS=$threads
    test_command test_border_values lh.surf mri_brain.mgz mri_smooth.mgz mri_aseg.mgz
done
