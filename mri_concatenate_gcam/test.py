#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

# copy LTA
rt.run('mri_concatenate_gcam im1_to_im2.lta copy.lta')
rt.diff('copy.lta', 'copy.ref.lta', ignore_comments=True)

# concatenate LTAs
rt.run('mri_concatenate_gcam im3_to_im1.lta im1_to_im2.lta concat.lta')
rt.diff('concat.lta', 'concat.ref.lta', ignore_comments=True)

# make sure we have mri_convert available before continuing
mri_convert = rt.findPath(rt.testdatadir, 'mri_convert/mri_convert')

# change the GCAM source image space
rt.run('mri_concatenate_gcam -s im1.mgz im2_to_im3.m3z gcam.m3z && ' +
        mri_convert + ' -at gcam.m3z im2_in_space1.mgz src.mgz')
rt.mridiff('src.mgz', 'src.ref.mgz')

# change the GCAM target image space
rt.run('mri_concatenate_gcam -t im1.mgz im2_to_im3.m3z gcam.m3z && ' +
        mri_convert + ' -at gcam.m3z im2.mgz dest.mgz')
rt.mridiff('dest.mgz', 'dest.ref.mgz')

# change the GCAM target image space to a larger space
rt.run('mri_concatenate_gcam --change-target im2.mgz im2_to_im1.m3z gcam.m3z')
rt.mridiff('gcam.m3z', 'larger.ref.m3z')

# concatenate LTA with GCAM
rt.run('mri_concatenate_gcam im1_to_im2.lta im2_to_im3.m3z gcam.m3z && ' +
        mri_convert + ' -at gcam.m3z im1.mgz lta_gcam.mgz')
rt.mridiff('lta_gcam.mgz', 'lta_gcam.ref.mgz')

# concatenate GCAM with LTA
rt.run('mri_concatenate_gcam im2_to_im3.m3z im3_to_im1.lta gcam.m3z && ' +
       mri_convert + ' -at gcam.m3z im2.mgz gcam_lta.mgz')
rt.mridiff('gcam_lta.mgz', 'gcam_lta.ref.mgz')

# concatenate LTA with GCAM with LTA
rt.run('mri_concatenate_gcam im1_to_im2.lta im2_to_im3.m3z im3_to_im1.lta gcam.m3z && ' +
       mri_convert + ' -at gcam.m3z im1.mgz lta_gcam_lta.mgz')
rt.mridiff('lta_gcam_lta.mgz', 'lta_gcam_lta.ref.mgz')

# concatenate two GCAMs
rt.run('mri_concatenate_gcam im2_to_im3.m3z im3_to_im2.m3z gcam.m3z && ' +
       mri_convert + ' -at gcam.m3z im2.mgz gcam_gcam.mgz')
rt.mridiff('gcam_gcam.mgz', 'gcam_gcam.ref.mgz')

# downsample GCAM - leave some capacity for change
rt.run('mri_concatenate_gcam -d im2_to_im3.m3z gcam.m3z && ' +
       mri_convert + ' -at gcam.m3z im2.mgz down.mgz')
rt.mridiff('down.mgz', 'down.ref.mgz', thresh=0.5)

rt.cleanup()

