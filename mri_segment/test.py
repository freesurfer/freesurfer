#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mri_segment mri/brainmask.mgz wm_segment_out.mgz -wlo 80')
rt.mridiff('wm_segment_out.mgz', 'wm_segment_ref.mgz')

rt.cleanup()
