#!/usr/bin/env python
import sys, os, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst
import platform

# platform has become distro in python 3
if sys.version_info[0] >= 3:
    import distro
else:
    import platform as distro

platform = distro.system()
os_params = ""
if (platform == "Linux"):
   os_params = distro.linux_distribution()
   # e.g., for CentOS 7 the returned tuple is , ('CentOS Linux', '7.4.1708', 'Core')
   os_list = os_params[0].split()
   os_name = os_list[0]
   os_rev = os_params[1]
   
rt = fst.RegressionTest()

rt.run('.//mri_surf2vol --o out.nii.gz --subject cvs_avg35 --so lh.white lh.thickness --so rh.white rh.thickness')

if (platform == "Darwin"):
   print ("Running diff for Mac OS")
   rt.diff('out.nii.gz', 'out_darwin_ref.nii.gz')
elif (platform == "Linux") and (os_name == "CentOS") and (os_rev[0] == "6"):
   print ("Running diff for Cent OS6")
   rt.diff('out.nii.gz', 'out_linux_CentOS6_ref.nii.gz')
elif (platform == "Linux") and (os_name == "CentOS") and (os_rev[0] == "7"):
   print ("Running diff for Cent OS7")
   rt.diff('out.nii.gz', 'out_linux_CentOS7_ref.nii.gz')
else:
   print ("Need more data to run on this platform:")
   if (platform == "Linux"): print ("%s\n") % str(os_params)
   sys.exit(1)

rt.cleanup()

