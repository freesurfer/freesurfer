#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mri_diff lh.ribbon.mgz rh.ribbon.mgz --ssd --count > diff.log', allow_failure=True)
rt.diff('diff.log', 'diff_ref.log')

rt.cleanup()
