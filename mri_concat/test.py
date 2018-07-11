#!/usr/bin/env python
import os
import sys, os.path as op
import os
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst
import glob

home = os.environ.get("HOME")
home_license = home + "/" + "license.txt"
if not os.environ.get("FS_LICENSE") and op.isfile(home_license): 
   os.environ["FS_LICENSE"] = home_license

rt = fst.RegressionTest()

rt.run('./mri_concat std.rh.*.mgh --o rhout.mgh')
rt.mridiff('rhout.mgh', 'rhout.ref.mgh')

rt.cleanup()

