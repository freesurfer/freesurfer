#!/usr/bin/env python
import sys, os.path as op, os
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

bdir = str(os.getcwd()) + '/testdata'
os.environ['SUBJECTS_DIR'] = bdir

grep_pattern = 'grep -v "^reading"'

cmd_lh = 'mris_surface_stats -mask ' + bdir + '/bert/label/lh.cortex.label -nsmooth 60 -surf_name ' + bdir + '/bert/surf/lh.white -src_type paint -out_name lh_std_60.mgh -absmean lh_absmean_60.mgh -mean lh_mean_60.mgh -absstd lh_absstd_60.mgh ' + bdir + '/bert/surf/lh.thickness 2>&1 | tee -a lh_surface_stats.out.raw'
cmd_lh_filter = 'cat lh_surface_stats.out.raw | ' + grep_pattern + ' > lh_surface_stats.out'
cmd_lh_all = cmd_lh + ' && ' + cmd_lh_filter

cmd_rh = 'mris_surface_stats -mask ' + bdir + '/bert/label/rh.cortex.label -nsmooth 60 -surf_name ' + bdir + '/bert/surf/rh.white -src_type paint -out_name rh_std_60.mgh -absmean rh_absmean_60.mgh -mean rh_mean_60.mgh -absstd rh_absstd_60.mgh ' + bdir + '/bert/surf/rh.thickness 2>&1 | tee -a rh_surface_stats.out.raw'
cmd_rh_filter = 'cat rh_surface_stats.out.raw | ' + grep_pattern + ' > rh_surface_stats.out'
cmd_rh_all = cmd_rh + ' && ' + cmd_rh_filter

cmd_all = cmd_lh_all + " && " + cmd_rh_all

rt.run(cmd_all)

rt.diff('lh_surface_stats.out', 'lh_surface_stats.out.ref')
rt.diff('rh_surface_stats.out', 'rh_surface_stats.out.ref')

rt.cleanup()

