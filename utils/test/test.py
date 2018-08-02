#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('testcolortab %s/FreeSurferColorLUT.txt' % rt.fs_home)
rt.run('test_c_nr_wrapper')
rt.run('extest')
rt.run('inftest')
rt.run('tiff_write_image')
rt.run('sc_test')

rt.cleanup()
