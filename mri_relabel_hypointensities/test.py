#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mri_relabel_hypointensities mri/aseg.presurf.mgz surf aseg_hypoint_new.mgz')
rt.mridiff('aseg_hypoint_new.mgz', 'aseg_hypoint_ref.mgz')

rt.cleanup()
