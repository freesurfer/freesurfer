#!/usr/bin/env python
import sys, os, os.path as op, platform, re
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

# set these to False for automated runs
# diff_only = True
diff_only = False
# manual_diff = True
manual_diff = False

# Darwin or Linux
platform = platform.system()
ref_file = "aseg.auto_noCCseg.mgz"
if (platform == "Darwin"):
   mac_os_rev=os.system("sw_vers -productVersion")
   # ... for older versions of Mac OS
   # Lion, Mountain Lion 
   if re.search(r'^10\.[7-8]',str(mac_os_rev)): ref_file = aseg.auto_noCCseg.macos.lion.mgz
   # Mavericks
   if re.search(r'^10\.9',str(mac_os_rev)): ref_file = aseg.auto_noCCseg.macos.mavericks.mgz
   # Yosemite, El Capitan
   if re.search(r'^10\.[10-11]',str(mac_os_rev)): ref_file = aseg.auto_noCCseg.macos.el_cap.mgz

init_thread_cnt = '8'
thread_cntlist = [init_thread_cnt]

n_arch = "norm.mgz"
t_arch = "talairach.m3z"
rb_arch = "../../distribution/average/RB_all_2016-05-10.vc700.gca"

rt = fst.RegressionTest()

for thread_cnt in thread_cntlist:  
    os.environ['OMP_NUM_THREADS'] = thread_cnt
    aseg_arch = "aseg.auto_noCCseg." + thread_cnt + "cpu.mgz"
    print "thread_cnt = %s" % (thread_cnt)
    test_cmd = 'mri_ca_label' + ' ' \
       + '-relabel_unlikely 9 .3' + ' ' \
       + '-prior 0.5' + ' ' \
       + '-align' + ' ' \
       + n_arch + ' ' + t_arch + ' ' + rb_arch + ' ' + aseg_arch
    if not diff_only:
       print "Will run test cmd: %s" % (test_cmd)
       rt.run(test_cmd)

    if manual_diff:
       # check with previous, "manual", diff, i.e., cd there and compare to ref file
       diff_cmd = '../../mri_diff/mri_diff' + ' ' + ref_file + ' ' + aseg_arch + ' ' + '--debug' + ' ' + '--diff' + ' ' + ref_file
       print "Will run diff cmd: %s" % (diff_cmd)
       os.chdir("testdata")
       os.system(diff_cmd)
    else:
       rt.mridiff(ref_file,aseg_arch)

# don't delete test results if want to run diff command again
if not diff_only: rt.cleanup()

