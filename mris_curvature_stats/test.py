#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst
import platform

platform = platform.system()
rt = fst.RegressionTest()

rt.run('mris_curvature_stats -h 10 -G cvs_avg35 rh rh.curv 2>&1 | tee -a rh.curv.out && \
mris_curvature_stats -h 10 -G cvs_avg35 lh lh.curv 2>&1 | tee -a lh.curv.out')

if (platform == "Linux"):
   rt.diff('rh.curv.out', 'rh.curv.linux.out.ref')
   rt.diff('lh.curv.out', 'lh.curv.linux.out.ref')
elif (platform == "Darwin"):
   rt.diff('rh.curv.out', 'rh.curv.darwin.out.ref')
   rt.diff('lh.curv.out', 'lh.curv.darwin.out.ref')

rt.cleanup()

