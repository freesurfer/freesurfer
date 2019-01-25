#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

test_command mri_surfcluster --in bert/surf/rh.w-g.pct.mgh --hemi rh --subject bert --thmin 22 --ocn rh.clusters.mgh --olab cluster
compare_vol rh.clusters.mgh rh.clusters_ref.mgh

for label in ${FSTEST_TESTDATA_DIR}/bert/label/cluster-*.label; do 
  compare_file $label ${label}.ref
done

