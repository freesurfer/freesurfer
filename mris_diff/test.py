#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mris_diff lh.sphere lh.sphere')
rt.run('mris_diff lh.sphere rh.sphere', expect_failure=True)

rt.cleanup()
