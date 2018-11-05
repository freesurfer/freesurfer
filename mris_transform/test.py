#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

# transform surface with LTA
rt.run('mris_transform lh.1.pial 1_to_2.lta lh.out.pial')
rt.surfdiff('lh.out.pial', 'lh.ref.pial')

# transform surface with inverse LTA
rt.run('mris_transform --is-inverse lh.1.pial 2_to_1.lta lh.out.pial')
rt.surfdiff('lh.out.pial', 'lh.ref.pial')

# transform surface with FSLMAT
rt.run(('mris_transform'
        ' --trx-src nu.1.mgz'
        ' --trx-dst nu.2.mgz'
        ' lh.1.pial 1_to_2.fslmat lh.out.pial'))
rt.surfdiff('lh.out.pial', 'lh.ref.pial')

# transform surface with FSLMAT
rt.run(('mris_transform'
        ' --trx-src nu.2.mgz'
        ' --trx-dst nu.1.mgz'
        ' --is-inverse'
        ' lh.1.pial 2_to_1.fslmat lh.out.pial'))
rt.surfdiff('lh.out.pial', 'lh.ref.pial')

# transform surface with GCAM
rt.run(('mris_transform'
        ' --trx-dst nu.2.mgz'
        ' lh.1.pial 1_to_2.m3z lh.out.pial'))
rt.surfdiff('lh.out.pial', 'lh.ref.pial')

# transform surface with inverse GCAM
rt.run(('mris_transform'
        ' --trx-dst nu.1.mgz'
        ' --is-inverse'
        ' lh.1.pial 2_to_1.m3z lh.out.pial'))
rt.surfdiff('lh.out.pial', 'lh.ref.pial')

rt.cleanup()

