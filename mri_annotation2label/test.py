#!/usr/bin/env python
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

rt.run('./mri_annotation2label --subject ./bert --sd ./. --hemi rh --outdir ./labels_from_test_run')

# cd'd to ./testdata
# edit files before diff'infg - ref and newly produced files from test have the same names - but in different subdirs
for filename in glob.iglob('./labels_from_test_run/*.label'):
     f_ref = "./labels_ref/" + op.basename(filename)
     f_new = "./labels_from_test_run/" + op.basename(filename)
     f_edit = f_new  + ".edit"
     # diff after first line
     cmd1 = "tail -n +2 " + f_new + " > " + f_edit
     # print ("Running: %s" % (cmd1))
     os.system(cmd1)

     cmd2 = "cp -p -f " + f_edit + " " + f_new
     # print ("Running: %s" % (cmd2))
     os.system(cmd2)

     cmd3 = "rm -f " + f_edit
     # print ("Running: %s" % (cmd3))
     os.system(cmd3)
     rt.diff(f_new,f_ref)

rt.cleanup()

