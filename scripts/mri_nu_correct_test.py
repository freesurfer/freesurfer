#!/usr/bin/env python
import sys, os, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mri_nu_correct.mni --i rawavg.mgz --o output.mgz --n 4 --no-uchar')
rt.mridiff('output.mgz', 'output.ref.mgz')

rt.cleanup()
