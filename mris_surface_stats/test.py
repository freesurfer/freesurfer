#!/usr/bin/env python
import sys, os, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mris_surface_stats -mask bert/label/lh.cortex.label -nsmooth 60 -surf_name bert/surf/lh.white bert/surf/lh.thickness &> out.log')
rt.diff('out.log', 'ref.log')

rt.cleanup()
