import os
import os.path as op
import sys
import argparse
from .util import run, rmdir, rmext
from .log import term, errorExit


# This class is built specifically for configuring regression tests on 
# freesurfer binaries within the main source tree. This assumes all test data is
# stored as a 'testdata.tar.gz' tarball in the directory where the test script
# lives. Commands should be tested using run(), and outputs should be compared
# against reference data by using the diff(), mridiff, or surfdiff() functions.
# Example:
#
#   rt = fst.RegressionTest()
#
#   rt.run('mri_normalize -mprage nu.mgz T1.mgz')
#   rt.mridiff('T1.mgz', 'T1.ref.mgz')
#
#   rt.run('mri_normalize -gentle orig.mgz gentle.mgz')
#   rt.mridiff('gentle.mgz', 'gentle.ref.mgz')
#
#   rt.cleanup()
#
# The initialization of a RegressionTest object will parse the command line args and
# look for the --regenerate flag, which will overwrite the reference data with any 
# produced outputs specified in the diff functions.

class RegressionTest:
  def __init__(self):
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--regenerate', action='store_true', help='regenerate the reference data')
    parser.add_argument('--keep-data', action='store_true', help='keep the testdata dir even after success')
    args = parser.parse_args()
    self.keepdata = args.keep_data
    # get directory of the test script
    self.testdir = os.getcwd()
    self.testdatadir = op.join(self.testdir, 'testdata')
    self.scriptdir = op.dirname(op.realpath(sys.argv[0]))
    self.testdatatar = op.join(self.scriptdir, 'testdata.tar.gz')
    # set up FREESURFER_HOME and default SUBJECTS_DIR
    fs_home = self.findPath(self.scriptdir, 'distribution')
    os.environ['FREESURFER_HOME'] = fs_home
    os.environ['SUBJECTS_DIR'] = self.testdatadir
    # set up martinos license if needed
    have_license = os.environ.get('FS_LICENSE') or op.exists(op.join(fs_home, '.license')) or op.exists(op.join(fs_home, 'license.txt'))
    if not have_license and op.exists('/space/freesurfer/.license'):
      os.environ["FS_LICENSE"] = '/space/freesurfer/.license'
    # if regenerating...
    self.regenerate = args.regenerate
    if self.regenerate:
      print('regenerating %s testdata' % op.basename(self.scriptdir))
      # make the tmp generation directory
      self.regeneration_dir = op.join(self.testdir, 'testdata_regeneration')
      rmdir(self.regeneration_dir)
      os.makedirs(self.regeneration_dir)
      os.chdir(self.regeneration_dir)
      # extract the original testdata to get replaced later
      self.runcmd('tar -xzvf ' + self.testdatatar)
      os.chdir(self.testdir)


  # recursively search upwards for a given target
  # this is used to locate the 'distribution' folder or the compiled 'mri_diff'
  # in the build directory
  def findPath(self, startpath, target):
    distpath = op.join(startpath, target)
    if not op.exists(distpath):
      parent = op.abspath(op.join(startpath, '..'))
      if parent == '/': errorExit("max recursion reached - could not locate '%s'" % target)
      return self.findPath(parent, target)
    return distpath


  # print and run a command, and exit if it fails
  def runcmd(self, cmd, fatal=True):
    print(term.bold + '>> ' + cmd + term.end)
    ret = run(cmd)
    if fatal and ret != 0:
      errorExit('"%s" failed' % cmd)
    return ret


  # a fatal cd
  def cd(self, dirname):
    if not op.isdir(dirname): errorExit('directory "%s" does not exist' % dirname)
    os.chdir(dirname)


  # run a test command (this is what should be called from the actual test script)
  def run(self, cmd, threads=None, allow_failure=False):
    self.cd(self.testdir)
    # extract the testdata
    rmdir(self.testdatadir)
    self.runcmd('tar -xzvf ' + self.testdatatar)
    self.cd(self.testdatadir)
    # set number of OMP threads if necessary
    if threads is not None:
      os.environ['OMP_NUM_THREADS'] = str(threads)
      print('setting OMP_NUM_THREADS = %s' % os.environ['OMP_NUM_THREADS'])
    # run the test command - assumes the command is located in the testdir (above testdatadir)
    if self.runcmd('../' + cmd, fatal=False) != 0 and not allow_failure:
      errorExit('test command "%s" failed' % cmd)


  # a simple diff using the standard unix diff command
  def diff(self, orig, ref):
    if self.regenerate:
      self._regen(orig, ref)
    else:
      self.cd(self.testdatadir)
      cmd = 'diff %s %s' % (orig, ref)
      if self.runcmd(cmd, fatal=False) != 0:
        errorExit('diff of %s and %s failed' % (orig, ref))


  # run a diff on two volumes (calls mri_diff, which must be already
  # built in the source directory)
  def mridiff(self, orig, ref, thresh=0.0, res_thresh=1e-6, geo_thresh=8e-6, flags=""):
    if self.regenerate:
      self._regen(orig, ref)
    else:
      self.cd(self.testdatadir)
      diffcmd = op.relpath(self.findPath(self.testdatadir, 'mri_diff/mri_diff'))
      cmd = '%s %s %s' % (diffcmd, orig, ref)
      origname = rmext(op.basename(orig))
      cmd += ' --debug --diff diff-%s.mgz --log diff-%s.log' % (origname, origname)
      cmd += ' --thresh %f --res-thresh %f --geo-thresh %f %s' % (thresh, res_thresh, geo_thresh, flags)
      if self.runcmd(cmd, fatal=False) != 0:
        errorExit('mri_diff of %s and %s failed' % (orig, ref))


  # run a diff on two surfs (calls mris_diff, which must be already
  # built in the source directory)
  def surfdiff(self, orig, ref, flags=""):
    if self.regenerate:
      self._regen(orig, ref)
    else:
      self.cd(self.testdatadir)
      diffcmd = op.relpath(self.findPath(self.testdatadir, 'mris_diff/mris_diff'))
      cmd = '%s %s %s --debug %s' % (diffcmd, orig, ref, flags)
      if self.runcmd(cmd, fatal=False) != 0:
        errorExit('mris_diff of %s and %s failed' % (orig, ref))
      self.cd(self.testdir)  # return to testdir


  # overwrite a reference file with the observed output
  # this will re-tar and overwrite the primary testdata.tar.gz
  def _regen(self, orig, ref):
    os.chdir(self.regeneration_dir)
    # replace the old with the new
    self.runcmd('mv -f %s %s' % (op.join(op.relpath(self.testdatadir), orig), op.join('testdata', ref)))
    self.cd(self.testdatadir)


  # delete the testdata directory and output any final messages
  # this should only get called when the whole test finishes successfully
  def cleanup(self):
    self.cd(self.testdir)
    if not self.keepdata: rmdir('testdata')
    # tar the regenerated data
    if self.regenerate:
      # tar up the testdata
      self.cd(self.regeneration_dir)
      self.runcmd('tar -czvf testdata.tar.gz testdata')
      # make sure the annex file is unlocked before replacing it
      os.chdir(self.scriptdir)
      self.runcmd('git annex unlock testdata.tar.gz')
      self.runcmd('mv -f %s .' % op.join(self.regeneration_dir, 'testdata.tar.gz'))
      rmdir(self.regeneration_dir)
      print('testdata has been regenerated')
      print('make sure to run "git annex add testdata.tar.gz" to rehash before committing')
      self.cd(self.testdir)
    else:
      print('%s[success]%s %s test has passed' % (term.green, term.end, op.basename(self.scriptdir)))
