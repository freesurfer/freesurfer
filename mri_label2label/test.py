#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mri_label2label --srclabel bert/label/lh.BA1_exvivo.label --srcsubject bert --trglabel cvs_avg35/label/lh.BA1.label --trgsubject cvs_avg35  --regmethod surface --hemi lh')
rt.diff('cvs_avg35/label/lh.BA1.label', 'label_ref/lh.BA1.label')

rt.cleanup()
