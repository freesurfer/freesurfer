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

# cmd = "FOO=/a/b/c date"
# ... results in ...
# >> ../FOO=/a/b/c date

# cmd = "SUBJECTS_DIR=" + subjects_dir + " mri_label2label --srclabel $SUBJECTS_DIR/bert/label/lh.BA1_exvivo.label  --srcsubject bert --trglabel $SUBJECTS_DIR/cvs_avg35/label/lh.BA1.label --trgsubject cvs_avg35  --regmethod surface --hemi lh"

cmd = "mri_label2label --srclabel $SUBJECTS_DIR/bert/label/lh.BA1_exvivo.label  --srcsubject bert --trglabel $SUBJECTS_DIR/cvs_avg35/label/lh.BA1.label --trgsubject cvs_avg35  --regmethod surface --hemi lh"

rt.run(cmd)

f_out = subjects_dir + "/cvs_avg35/label/lh.BA1.label"
f_ref = subjects_dir + "/label_ref/lh.BA1.label" 
rt.diff(f_out,f_ref)

rt.cleanup()

