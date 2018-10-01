#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mri_binarize --i aseg.mgz --o mri_binarize.mgz --match 17')
rt.mridiff('mri_binarize.mgz', 'mri_binarize.ref.mgz')

rt.cleanup()
