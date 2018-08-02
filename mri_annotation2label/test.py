#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst
import glob

rt = fst.RegressionTest()

rt.run('mri_annotation2label --subject bert --sd . --hemi rh --outdir labels')

for filename in glob.iglob('labels/*.label'):
  rt.diff(filename, op.join('labels_ref', op.basename(filename)))

rt.cleanup()
