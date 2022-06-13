#!/usr/bin/env bash
source "$(dirname $0)/../test.sh"

# conform
test_command mri_convert rawavg.mgz orig.mgz --conform
compare_vol orig.mgz orig.ref.mgz

# dicom
test_command mri_convert dcm/261000-10-60.dcm dicom.mgz
compare_vol dicom.mgz freesurfer.mgz

# nifti
test_command mri_convert nifti.nii nifti.mgz
compare_vol nifti.mgz freesurfer.mgz --notallow-acq --geo-thresh 0.000008

# analyze
test_command mri_convert analyze.img analyze.mgz
compare_vol analyze.mgz freesurfer.mgz --notallow-acq --geo-thresh 0.000008

# mri_make_uchar
test_command mri_make_uchar nu.mgz talairach.xfm nuuc.mgz
compare_vol nuuc.mgz nuuc.ref.mgz

# apply downsampled morph with atlas geometry that has odd number of slices
# (i.e. gcam->depth * gcam->spacing = gcam->atlas.depth - 1)
test_command mri_convert -at odd.m3z orig.mgz morphed.mgz
compare_vol morphed.mgz odd.ref.mgz

# standard mosaic'd DICOM
test_command mri_convert ep2d.mosaic.dcm ep2d.mosaic.mgz
compare_vol ep2d.mosaic.mgz ep2d.mosaic.ref.mgz

# non-mosaic DICOM with incomplete ASCII header
test_command mri_convert vnav.non-mosaic.dcm vnav.non-mosaic.mgz
compare_vol vnav.non-mosaic.mgz vnav.non-mosaic.ref.mgz

# DICOM with identical geometry - but mosaic'd
test_command mri_convert --mosaic-fix-noascii vnav.mosaic.dcm vnav.mosaic.mgz
compare_vol vnav.mosaic.mgz vnav.mosaic.ref.mgz

# -dcm2niix DICOM conversion - kSliceOrientMosaicNegativeDeterminant (the data is taken from dcm2niix qa data)
test_command mri_convert -dcm2niix dtitest_Siemens_ccbbi/DTI_sag_002_001_00001.dcm dtinegativemosaic.mgz
compare_vol dtinegativemosaic.mgz dtitest_Siemens_ccbbi/ref/dtinegativemosaic.ref.mgz
compare_file dtinegativemosaic.bvals dtitest_Siemens_ccbbi/ref/dtinegativemosaic.ref.bvals
compare_file dtinegativemosaic.voxel_space.bvecs dtitest_Siemens_ccbbi/ref/dtinegativemosaic.ref.voxel_space.bvecs

# -dcm2niix DICOM conversion - sliceDir < 0 
test_command mri_convert -dcm2niix enhanced28/34983475 enhanced28.mgz
compare_vol enhanced28.mgz enhanced28/ref/enhanced28.ref.mgz
compare_file enhanced28.bvals enhanced28/ref/enhanced28.ref.bvals
compare_file enhanced28.voxel_space.bvecs enhanced28/ref/enhanced28.ref.voxel_space.bvecs

# -dcm2niix DICOM conversion - dti default case
test_command mri_convert -dcm2niix enhanced35/34984333 enhanced35.mgz
compare_vol enhanced35.mgz enhanced35/ref/enhanced35.ref.mgz
compare_file enhanced35.bvals enhanced35/ref/enhanced35.ref.bvals
compare_file enhanced35.voxel_space.bvecs enhanced35/ref/enhanced35.ref.voxel_space.bvecs
