#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst
import glob

rt = fst.RegressionTest()

rt.run('mri_surfcluster --in bert/surf/rh.w-g.pct.mgh --hemi rh --subject bert --thmin 22 --ocn rh.clusters.mgh --olab cluster')
rt.mridiff('rh.clusters.mgh', 'rh.clusters_ref.mgh')
for label in glob.iglob('bert/label/cluster-*.label'):
  rt.diff(label, label + '.ref')

rt.cleanup()
