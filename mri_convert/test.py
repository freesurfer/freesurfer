#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

# conform
rt.run('mri_convert rawavg.mgz orig.mgz --conform')
rt.mridiff('orig.mgz', 'orig.ref.mgz')

# dicom
rt.run('mri_convert dcm/261000-10-60.dcm dicom.mgz')
rt.mridiff('dicom.mgz', 'freesurfer.mgz')

# nifti
rt.run('mri_convert nifti.nii nifti.mgz')
rt.mridiff('nifti.mgz', 'freesurfer.mgz', flags='--notallow-acq')

# analyze
rt.run('mri_convert analyze.img analyze.mgz')
rt.mridiff('analyze.mgz', 'freesurfer.mgz', flags='--notallow-acq')

# mri_make_uchar
rt.run('mri_make_uchar nu.mgz talairach.xfm nuuc.mgz')
rt.mridiff('nuuc.mgz', 'nuuc.ref.mgz')

# apply downsampled morph with atlas geometry that has odd number of slices
# (i.e. gcam->depth * gcam->spacing = gcam->atlas.depth - 1)
rt.run('mri_convert -at odd.m3z orig.mgz morphed.mgz')
rt.mridiff('morphed.mgz', 'odd.ref.mgz')

rt.cleanup()
