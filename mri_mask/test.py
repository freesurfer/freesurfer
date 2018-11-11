#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

# apply mask
rt.run('mri_mask nu.1.mgz brainmask.1.mgz mask.mgz')
rt.mridiff('mask.mgz', 'mask.ref.mgz')

# transform mask using LTA
rt.run('mri_mask -xform 2_to_1.lta nu.1.mgz brainmask.2.mgz mask.lta.mgz')
rt.mridiff('mask.lta.mgz', 'mask.lta.ref.mgz')

# transform mask using inverse LTA
rt.run('mri_mask -xform 1_to_2.lta nu.1.mgz brainmask.2.mgz mask.lta.inv.mgz')
rt.mridiff('mask.lta.inv.mgz', 'mask.lta.inv.ref.mgz')

# transform mask with FSLMAT
rt.run(('mri_mask'
        ' -xform 2_to_1.fslmat'
        ' -lta_src nu.2.mgz'
        ' -lta_dst nu.1.mgz'
        ' nu.1.mgz brainmask.2.mgz mask.fsl.mgz'))
rt.mridiff('mask.fsl.mgz', 'mask.fsl.ref.mgz')

# transform mask with GCAM
rt.run('mri_mask -xform 2_to_1.m3z nu.1.mgz brainmask.2.mgz mask.gcam.mgz')
rt.mridiff('mask.gcam.mgz', 'mask.gcam.ref.mgz')

# transfer edits
rt.run(('mri_mask'
        ' -transfer 255'
        ' -keep_mask_deletion_edits'
        ' nu.2.mgz wm.2.mgz edits.mgz'))
rt.mridiff('edits.mgz', 'edits.ref.mgz')

# transfer edits using identity.nofile
rt.run(('mri_mask'
        ' -xform identity.nofile'
        ' -transfer 255'
        ' -keep_mask_deletion_edits'
        ' nu.2.mgz wm.2.mgz edits.mgz'))
rt.mridiff('edits.mgz', 'edits.ref.mgz')

# transfer edits using LTA
rt.run(('mri_mask'
        ' -xform 2_to_1.lta'
        ' -transfer 255'
        ' -keep_mask_deletion_edits'
        ' nu.1.mgz wm.2.mgz edits.lta.mgz'))
rt.mridiff('edits.lta.mgz', 'edits.lta.ref.mgz')

# transfer edits using GCAM
rt.run(('mri_mask'
        ' -xform 2_to_1.m3z'
        ' -transfer 255'
        ' -keep_mask_deletion_edits'
        ' nu.1.mgz wm.2.mgz edits.gcam.mgz'))
rt.mridiff('edits.gcam.mgz', 'edits.gcam.ref.mgz')

# transfer 255 only using inverse GCAM
rt.run(('mri_mask'
        ' -xform 1_to_2.m3z'
        ' -transfer 255'
        ' nu.1.mgz wm.2.mgz edits.gcam.inv.mgz'))
rt.mridiff('edits.gcam.inv.mgz', 'edits.gcam.inv.ref.mgz')

rt.cleanup()

