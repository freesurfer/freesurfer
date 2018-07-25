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

# should be <abspath>/testdata
subjects_dir = os.getcwd() + "/testdata"
os.environ["SUBJECTS_DIR"] = subjects_dir

rt.run('./mri_fuse_segmentations -a aseg.mgz -c aseg.auto_noCCseg.mgz -n norm.mgz bert fused.mgz')

# sys.exit(0)

rt.mridiff("fused.mgz", "fused.ref.mgz")

rt.cleanup()

