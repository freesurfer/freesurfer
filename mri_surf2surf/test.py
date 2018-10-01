#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mri_surf2surf --hemi lh --srcsubject bert --srcsurfval thickness --src_type curv --trgsubject fsaverage --trg_type curv --trgsurfval bert/surf/lh.thickness')
rt.mridiff('bert/surf/lh.thickness', 'bert/surf/lh.ref_thickness')

rt.cleanup()
