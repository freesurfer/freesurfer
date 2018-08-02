#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mri_matrix_multiply -v -im bert/mri/transforms/talairach.xfm -im $SUBJECTS_DIR/cvs_avg35/mri/transforms/talairach.xfm -v -om mmult.xfm')
rt.diff('mmult.xfm', 'mmult.xfm.ref')

rt.cleanup()
