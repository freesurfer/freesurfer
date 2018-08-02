#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mri_pretess -w mri/wm.mgz wm mri/norm.mgz wm_new.mgz')
rt.mridiff('wm_new.mgz', 'wm_ref.mgz')

rt.cleanup()
