#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mri_fwhm --i HelixTensors.nii.gz --nframesmin 9 --auto-mask .2 --dat fwhm.dat')
rt.diff('fwhm.dat', 'fwhm_ref.dat')

rt.cleanup()
