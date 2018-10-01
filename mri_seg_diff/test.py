#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mri_seg_diff --debug --seg1 mri/aseg.auto.mgz --seg2 mri/aseg.mgz --diff aseg_diff_out.mgz')
rt.mridiff('aseg_diff_out.mgz', 'aseg_diff_ref.mgz')

rt.cleanup()
