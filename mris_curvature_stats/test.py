#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mris_curvature_stats -h 10 -G cvs_avg35 rh rh.curv 2>&1 | tee -a rh.curv.out && \
mris_curvature_stats -h 10 -G cvs_avg35 lh lh.curv 2>&1 | tee -a lh.curv.out')

rt.diff('rh.curv.out', 'rh.curv.out.ref')
rt.diff('lh.curv.out', 'lh.curv.out.ref')
rt.cleanup()

