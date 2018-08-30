#!/usr/bin/env python
import sys, os, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst
import platform

platform = platform.system()
rt = fst.RegressionTest()

rt.run('.//mri_surf2vol --o out.nii.gz --subject cvs_avg35 --so lh.white lh.thickness --so rh.white rh.thickness')

if (platform == "Linux"):
   rt.diff('out.nii.gz', 'out_linux_ref.nii.gz')
elif (platform == "Darwin"):
   rt.diff('out.nii.gz', 'out_darwin_ref.nii.gz')

rt.cleanup()
