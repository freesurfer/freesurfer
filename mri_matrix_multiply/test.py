#!/usr/bin/env python
import os
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

# should be <abspath>/testdata
subjects_dir = os.getcwd() + "/testdata"
os.environ["SUBJECTS_DIR"] = subjects_dir

home = os.environ.get("HOME")
home_license = home + "/" + "license.txt"
if not os.environ.get("FS_LICENSE") and op.isfile(home_license): 
   os.environ["FS_LICENSE"] = home_license

cmd = "mri_matrix_multiply -v -im $SUBJECTS_DIR/bert/mri/transforms/talairach.xfm -im $SUBJECTS_DIR/cvs_avg35/mri/transforms/talairach.xfm -v -om mmult.xfm"

rt.run(cmd)

f_out = subjects_dir + "/" + "mmult.xfm"
f_ref = subjects_dir + "/" + "mmult.xfm.ref" 
rt.diff(f_out,f_ref)

rt.cleanup()

