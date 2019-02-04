#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command optseq2 \
    --ntp 180 \
    --tr 2 \
    --psdwin 0 24 1 \
    --ev Neutral-Short 3 24 \
    --ev Neutral-Long  3 24 \
    --ev Fearful-Short 3 24 \
    --ev Fearful-Long  3 24 \
    --polyfit 2 \
    --tnullmax 10 \
    --focb 100 \
    --nsearch 100 \
    --nkeep 4 \
    --o emot \
    --seed 1234

# remove lines from output files that contain non-essential stuff
function filter {
    grep -v "optseq2.c,v" $1 > tmp && mv -f tmp $1
    grep -v "hours" $1 > tmp && mv -f tmp $1
    grep -v "iterations per second" $1 > tmp && mv -f tmp $1
}

# compare output
for f in emot-001.par emot-002.par emot-003.par emot-004.par emot.sum; do
    filter ${f}
    compare_file ${f} expected/${f}
done
