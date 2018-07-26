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

rt.run('./mri_fwhm --i HelixTensors.nii.gz --nframesmin 9 --auto-mask .2 --sum out.fwhm.sum')

# remove timestamp line whose date will always differ
f_new = "out.fwhm.sum"
f_edit = f_new  + ".edit"
cmd1 = 'cat ' + f_new + ' | grep -v ' + '"' + '^timestamp ' + '"' + ' > ' + f_edit
print ("Running: %s" % (cmd1))
os.system(cmd1)

rt.diff(f_edit, 'out.fwhm.sum.ref')

rt.cleanup()

