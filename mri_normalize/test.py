#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

# mprage
rt.run('mri_normalize -mprage nu.mgz T1.mgz')
rt.mridiff('T1.mgz', 'T1.ref.mgz')

# aseg
rt.run('mri_normalize -aseg aseg.presurf.mgz norm.mgz brain.mgz')
rt.mridiff('brain.mgz', 'brain.ref.mgz')

# gentle
rt.run('mri_normalize -gentle nu.mgz gentle.mgz')
rt.mridiff('gentle.mgz', 'gentle.ref.mgz')

rt.cleanup()
