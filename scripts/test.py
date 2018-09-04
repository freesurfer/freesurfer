#!/usr/bin/env python
import sys, os, os.path as op
import platform
sys.path.append(op.join(op.dirname(sys.argv[0]), '../python'))

home = os.environ['FREESURFER_HOME']
freesurfer_bin = home + '/bin'
freesurfer_mnibin = home  + '/mni/bin'
my_path = os.environ['PATH']
new_path = my_path + ':' + freesurfer_bin + ':' + freesurfer_mnibin
os.environ['PATH'] = new_path

platform = platform.system()
if (platform == "Darwin"):
   os.environ['PERL5LIB'] = home + '/mni/lib/../Library/Perl/Updates/5.12.3'
   os.environ['MNI_PERL5LIB'] = home + '/mni/lib/../Library/Perl/Updates/5.12.3'
elif (platform == "Linux"):
   os.environ['PERL5LIB'] = '/usr/pubsw/packages/mni/1.5/lib/perl5/5.8.5'
   os.environ['MNI_PERL5LIB'] = '/usr/pubsw/packages/mni/1.5/lib/perl5/5.8.5'

import freesurfer.test as fst

rt = fst.RegressionTest()

# mri_nu_correct.mni
rt.run('mri_nu_correct.mni --i rawavg.mgz --o output.mgz --n 4 --no-uchar')
rt.mridiff('output.mgz', 'output.ref.mgz')

rt.cleanup()
