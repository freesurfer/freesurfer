#!/usr/bin/env python
import os
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

home = os.environ.get("HOME")
home_license = home + "/" + "license.txt"
if not os.environ.get("FS_LICENSE") and op.isfile(home_license): 
   os.environ["FS_LICENSE"] = home_license

rt = fst.RegressionTest()

# mri_binarize
rt.run('mri_binarize --i aseg.mgz --o mri_binarize.mgz --match 17')
rt.mridiff('mri_binarize.mgz', 'mri_binarize.ref.mgz')

rt.cleanup()

