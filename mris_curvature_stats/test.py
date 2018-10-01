#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mris_curvature_stats -o raw -h 10 -G cvs_avg35 rh rh.curv')

rt.diff('raw', 'ref')
for stat in ('BE', 'C', 'FI', 'H', 'K1', 'K2', 'K', 'S', 'SI'):
  rt.diff('raw.rh.smoothwm.%s.hist' % stat, 'ref.rh.smoothwm.%s.hist' % stat)

rt.cleanup()
