#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mri_warp_convert --inm3z talairach.m3z --outitk out.nii.gz')
rt.mridiff('out.nii.gz', 'ref.nii.gz')

rt.cleanup()
