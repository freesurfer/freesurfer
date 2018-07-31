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

# remove timestamp, hostname, user, sysname, Id:, FREESURFER_HOME, cwd  which can differ depending upon where run
f_new = "out.fwhm.sum"
f_edit = f_new  + ".edit"
cmd1 = 'cat ' + f_new + ' | grep -v ' + '"' + '^timestamp ' + '"' + ' | grep -v ' + '"' + '^hostname ' + '"' + ' | grep -v ' + '"' + '^user ' + '"' + ' | grep -v ' + '"' + '^sysname ' + '"' + ' | grep -v ' + '"' + 'Id: ' + '"' +  ' | grep -v ' + '"' + '^FREESURFER_HOME ' + '"' + ' | grep -v ' + '"' + '^cwd ' + '"'  +  ' > ' + f_edit
print ("Running: %s" % (cmd1))
os.system(cmd1)

rt.diff(f_edit, 'out.fwhm.sum.ref')

rt.cleanup()

