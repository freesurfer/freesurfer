#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

# apply mask
rt.run('mri_mask nu.1.mgz brainmask.1.mgz mask.mgz')
rt.mridiff('mask.mgz', 'mask.ref.mgz')

# transform mask using LTA
rt.run('mri_mask -xform 2_to_1.lta nu.1.mgz brainmask.2.mgz lta.mgz')
rt.mridiff('lta.mgz', 'lta.ref.mgz')

# transform mask using inverse LTA
rt.run('mri_mask -xform 1_to_2.lta nu.1.mgz brainmask.2.mgz inv.mgz')
rt.mridiff('inv.mgz', 'inv.ref.mgz')

# transform mask with FSLMAT
rt.run(('mri_mask'
        ' -xform 2_to_1.fslmat'
        ' -lta_src nu.2.mgz'
        ' -lta_dst nu.1.mgz'
        ' nu.1.mgz brainmask.2.mgz fslmat.mgz'))
rt.mridiff('fslmat.mgz', 'fslmat.ref.mgz')

# transform mask with GCAM
rt.run('mri_mask -xform 2_to_1.m3z nu.1.mgz brainmask.2.mgz gcam.mgz')
rt.mridiff('gcam.mgz', 'gcam.ref.mgz')

# only transfer WM edits (255) and deletions (1)
rt.run(('mri_mask'
        ' -transfer 255'
        ' -keep_mask_deletion_edits'
        ' nu.2.mgz wm.2.mgz edits.mgz'))
rt.mridiff('edits.mgz', 'edits.ref.mgz')

rt.cleanup()

