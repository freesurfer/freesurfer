#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

# relabel corpus callosum
rt.run(('mri_fuse_segmentations'
        ' --aseg aseg.1.mgz'
        ' --nocc aseg.auto_noCCseg.1.mgz'
        ' --norm norm.1.mgz'
        ' norm.1.mgz cc.mgz'));
rt.mridiff('cc.mgz', 'cc.ref.mgz')

# fuse segmentations using LTAs
rt.run(('mri_fuse_segmentations'
        ' -a aseg.{1,2,3}.mgz'
        ' -c aseg.auto_noCCseg.{1,2,3}.mgz'
        ' -n norm.{1,2,3}.mgz'
        ' -t im{1,2,3}_to_base.lta'
        ' norm.base.mgz lta.mgz'));
rt.mridiff('lta.mgz', 'lta.ref.mgz')

# use different sigma
rt.run(('mri_fuse_segmentations'
        ' -a aseg.{1,2,3}.mgz'
        ' -c aseg.auto_noCCseg.{1,2,3}.mgz'
        ' -n norm.{1,2,3}.mgz'
        ' -t im{1,2,3}_to_base.lta'
        ' -s 4.0'
        ' norm.base.mgz sigma.mgz'));
rt.mridiff('sigma.mgz', 'sigma.ref.mgz')

# fuse segmentations using GCAMs and identity transform
rt.run(('mri_fuse_segmentations'
        ' -a aseg.{1,2,3}.mgz'
        ' -c aseg.auto_noCCseg.{1,2,3}.mgz'
        ' -n norm.{1,2,3}.mgz'
        ' -t im1_to_im2.m3z identity.nofile im3_to_im2.m3z'
        ' norm.2.mgz gcam.mgz'));
rt.mridiff('gcam.mgz', 'gcam.ref.mgz')

rt.cleanup()

