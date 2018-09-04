#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mri_fuse_segmentations -a aseg.mgz -c aseg.auto_noCCseg.mgz -n norm.mgz bert fused.mgz')
rt.mridiff("fused.mgz", "fused.ref.mgz")

rt.cleanup()

