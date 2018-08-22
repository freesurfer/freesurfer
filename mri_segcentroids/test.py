#!/usr/bin/env python
import sys, os, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('./mri_segcentroids --ctab-default --i aseg.mgz --o mri_segcentroids.raw')
# Not sure why -proper only sometimes appended to "Left => "Left-Thalamus"
# and Right => "Right-Thalmus" text entries, however numerical data is always the same
os.system("cat mri_segcentroids.raw | sed 's;-proper;;g' > mri_segcentroids.out")
rt.diff('mri_segcentroids.out', 'mri_segcentroids_ref.out')

rt.cleanup()
