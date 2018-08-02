#!/usr/bin/env python
import sys, os, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()
os.environ['FSDEV_TEST_DATA'] = rt.testdatadir
rt.run('test_fsgdf')
rt.cleanup()
