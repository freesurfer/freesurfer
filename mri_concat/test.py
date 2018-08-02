#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('./mri_concat std.rh.*.mgh --o rhout.mgh')
rt.mridiff('rhout.mgh', 'rhout.ref.mgh')

rt.cleanup()
