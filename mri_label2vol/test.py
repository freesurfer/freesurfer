#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mri_label2vol --label label/lh.cortex.label --regheader mri/rawavg.mgz --temp mri/aparc+aseg.mgz --o l2v.nii')
rt.mridiff('l2v.nii', 'l2v_ref.nii')

rt.cleanup()
