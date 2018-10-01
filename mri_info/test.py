#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst
import os, glob, re

rt = fst.RegressionTest()

string = ""
for filename in glob.iglob('./testdata/input/*.mgz'):
     f_in = re.sub('testdata/','',filename)
     f_new = "./output_from_test_run/" + op.basename(filename) + ".txt"
     cmd = "mri_info " + f_in + " > " + f_new
     if not string:
        string = cmd
     else:
        # only first command will have ../ appended since cd'd to ./testdata
        string = cmd + ' && ' + '../' + string

mri_info_abspath = os.path.abspath("../mri_info") + "/mri_info"
# prevent susbequent ../ with abspath to cmd
cmd_string = string.replace('../mri_info',mri_info_abspath)
print ("cmd_string = %s" % (cmd_string))

rt.run(cmd_string)

for filename in glob.iglob('./output_from_test_run/*.txt'):
     f_new = filename
     f_edit = filename + ".edit"
     f_ref = "./output_ref/" + op.basename(filename) + ".ref"

     cmd1 = "tail -n +2 " + f_new + " > " + f_edit
     print ("Running: %s" % (cmd1))
     os.system(cmd1)

     cmd2 = "cp -p -f " + f_edit + " " + f_new
     # print ("Running: %s" % (cmd2))
     os.system(cmd2)

     cmd3 = "rm -f " + f_edit
     # print ("Running: %s" % (cmd3))
     os.system(cmd3)
     rt.diff(f_new,f_ref)

# issue with shell and consecutive runs ?
# rt.cleanup()
