#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_glmfit \
    --seed 1234 \
    --y lh.gender_age.thickness.10.mgh \
    --fsgd gender_age.txt doss \
    --no-cortex \
    --glmdir lh.gender_age.glmdir \
    --surf average lh \
    --C age.mat

actual="${FSTEST_TESTDATA_DIR}/lh.gender_age.glmdir"
expected="${FSTEST_TESTDATA_DIR}/expected.lh.gender_age.glmdir"

compare_file ${actual}/Xg.dat ${expected}/Xg.dat
compare_file ${actual}/age/C.dat ${expected}/age/C.dat

for f in sar1.mgh beta.mgh mask.mgh rstd.mgh rvar.mgh; do
    compare_vol ${actual}/${f} ${expected}/${f} --thresh 0.00007
done

for f in F.mgh gamma.mgh sig.mgh; do
    compare_vol ${actual}/age/${f} ${expected}/age/${f} --thresh 0.008
done
