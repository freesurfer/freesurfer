#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

rt.run('mris_divide_parcellation bert lh aparc.annot splittable.txt split-sfg+pc 2>&1 | tee -a divide_parcel.out.raw && cat divide_parcel.out.raw | grep -v "supposed to be reproducible" > divide_parcel.out') 

rt.diff('divide_parcel.out', 'divide_parcel.out.ref')

rt.cleanup()

