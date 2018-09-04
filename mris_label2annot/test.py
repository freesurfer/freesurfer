#!/usr/bin/env python
import sys, os.path as op, os
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

bdir = str(os.getcwd()) + '/testdata'

cmd = 'mris_label2annot --sd ' + bdir + ' --s bert --h lh --ctab  bert/label/aparc.annot.ctab --a myaparc --l bert/label/lh.cortex.label --nhits nhits.mgh 2>&1 | tee -a label2annot.out'

rt.run(cmd)

rt.diff('bert/label/lh.myaparc.annot', 'bert/label/lh.myaparc.annot.ref')

rt.cleanup()

