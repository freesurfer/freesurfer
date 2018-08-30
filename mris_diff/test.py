#!/usr/bin/env python
import sys, os.path as op
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))
import freesurfer.test as fst

rt = fst.RegressionTest()

grep_pattern = 'grep -v "^timestamp" | grep -v "^sysname" | grep -v "^hostname" | grep -v "^machine" | grep -v "^user"'

# diff same surface
cmd_1 = 'mris_diff --debug --s1 cvs_avg35 ./cvs_avg35/surf/lh.sphere ./cvs_avg35/surf/lh.sphere 2>&1 | tee -a diff_1_subj_left_left.out.raw'
cmd_1_filter = 'cat diff_1_subj_left_left.out.raw | ' + grep_pattern + ' > diff_1_subj_left_left.out'
cmd_1_all = cmd_1 + ' && ' + cmd_1_filter

# diff left and right surface for same subject
cmd_2 = 'mris_diff --debug --s1 cvs_avg35 ./cvs_avg35/surf/lh.sphere ./cvs_avg35/surf/rh.sphere 2>&1 | tee -a diff_1_subj_left_right.out.raw'
cmd_2_filter = 'cat diff_1_subj_left_right.out.raw | ' + grep_pattern + ' > diff_1_subj_left_right.out'
cmd_2_all = cmd_2 + ' && ' + cmd_2_filter

# diff same surface for different subjects
cmd_3 = 'mris_diff --debug --s1 cvs_avg35 ./cvs_avg35/surf/lh.sphere --s2 bert ./bert/surf/lh.sphere 2>&1 | tee -a diff_2_subj_left_left.out.raw'
cmd_3_filter = 'cat diff_2_subj_left_left.out.raw | ' + grep_pattern + ' > diff_2_subj_left_left.out'
cmd_3_all = cmd_3 + ' && ' + cmd_3_filter

rt.run(cmd_1_all + '&&' + cmd_2_all + ' && ' + cmd_3_all)

# diff

rt.diff('diff_1_subj_left_left.out', 'diff_1_subj_left_left.out.ref')
rt.diff('diff_1_subj_left_right.out', 'diff_1_subj_left_right.out.ref')
rt.diff('diff_2_subj_left_left.out', 'diff_2_subj_left_left.out.ref')

rt.cleanup()

