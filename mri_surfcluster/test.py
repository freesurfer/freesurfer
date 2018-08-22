#!/usr/bin/env python
import sys, os.path as op
import os
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

os.environ['SUBJECTS_DIR'] = os.getcwd() + "/testdata"

rt.run('\
mri_surfcluster --in ./bert/surf/rh.w-g.pct.mgh --hemi rh --subject bert --thmin 0 2>&1 | tee -a mri_surfcluster.rh.out.raw && \
mri_surfcluster --in ./bert/surf/lh.w-g.pct.mgh --hemi lh --subject bert --thmin 0 2>&1 | tee -a mri_surfcluster.lh.out.raw && \
cat mri_surfcluster.rh.out.raw | \
grep -v "^[Vv]ersion" | grep -v "^[Rr]eading" | grep -v "^[Ss]ubjectsdir" | grep -v "^\/.*talairach.xfm"$ > mri_surfcluster.rh.out && \
cat mri_surfcluster.lh.out.raw | \
grep -v "^[Vv]ersion" | grep -v "^[Rr]eading" | grep -v "^[Ss]ubjectsdir" | grep -v "^\/.*talairach.xfm"$ > mri_surfcluster.lh.out \
')

rt.diff('mri_surfcluster.rh.out', 'mri_surfcluster.rh.out.ref')
rt.diff('mri_surfcluster.lh.out', 'mri_surfcluster.lh.out.ref')

rt.cleanup()

