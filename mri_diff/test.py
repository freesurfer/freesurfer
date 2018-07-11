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

f_ref = "ref_diff.txt"
f_new = "test_run_diff.txt"
# force exit 0 since testing diff result (returns non-zero status)
cmd = "mri_diff lh.ribbon.mgz rh.ribbon.mgz --log " + f_new + " ; exit 0" 
rt.run(cmd)

# edit output prior to diff, diff after 12th line
f_edit = f_new  + ".edit"
cmd1 = "tail -n +12 " + f_new + " > " + f_edit
print "Running: %s" % (cmd1)
os.system(cmd1)

cmd2 = "cp -p -f " + f_edit + " " + f_new
print "Running: %s" % (cmd2)
os.system(cmd2)

cmd3 = "rm -f " + f_edit
print "Running: %s" % (cmd3)
os.system(cmd3)

rt.diff(f_new,f_ref)
rt.cleanup()

