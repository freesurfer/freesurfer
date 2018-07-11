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

cmd = "mri_label2vol --label $SUBJECTS_DIR/label/lh.cortex.label --regheader $SUBJECTS_DIR/mri/rawavg.mgz --temp $SUBJECTS_DIR/mri/aparc+aseg.mgz --o $SUBJECTS_DIR/l2v.nii"

rt.run(cmd)

f_out = subjects_dir + "/" + "l2v.nii"
f_ref = subjects_dir + "/" + "l2v.nii.ref" 
rt.diff(f_out,f_ref)

rt.cleanup()

