import os
import os.path as op
import sys
import argparse
from . import run, rmdir, rmext, term, errorExit


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
    # get directory of the test script
    self.scriptdir = op.dirname(op.realpath(sys.argv[0]))
    self.testdatadir = op.join(self.scriptdir, 'testdata')
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--regenerate", action='store_true', help='regenerate the reference data')
    args = parser.parse_args()
    # setup FREESURFER_HOME
    os.environ["FREESURFER_HOME"] = self.findPath(self.scriptdir, 'distribution')
    # regenerate
    self.regenerate = args.regenerate
    if self.regenerate: print('regenerating %s testdata' % op.basename(self.scriptdir))


  # recursively search upwards for a given target
  # this is used to locate the 'distribution' folder or the compiled 'mri_diff'
  # in the main source directory
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
  # TODO: calls to this should probably just be replaced with self.runcmd('cd')
  def cd(self, dirname):
    if not op.isdir(dirname): errorExit('directory "%s" does not exist' % dirname)
    os.chdir(dirname)


  # run a test command (this is what should be called from the actual test script)
  def run(self, cmd):
    self.cd(self.scriptdir)
    # extract the testdata
    rmdir('testdata')
    self.runcmd('tar -xzvf testdata.tar.gz')
    self.cd('testdata')
    # run the test command (assumes the command is located above the testdata dir)
    if self.runcmd('../' + cmd, fatal=False) != 0:
      errorExit('test command "%s" failed' % cmd)


  # a simple diff using the standard unix diff command
  def diff(self, orig, ref):
    if self.regenerate:
      self.regen(orig, ref)
    else:
      self.cd(self.testdatadir)
      cmd = 'diff %s %s' % (orig, ref)
      if self.runcmd(cmd, fatal=False) != 0:
        errorExit('diff of %s and %s failed' % (orig, ref))


  # run a diff on two volumes (calls mri_diff, which must be already
  # built in the source directory)
  def mridiff(self, orig, ref, thresh=0.0, res_thresh=1e-6, geo_thresh=8e-6):
    if self.regenerate:
      self.regen(orig, ref)
    else:
      self.cd(self.testdatadir)
      diffcmd = op.relpath(self.findPath(self.testdatadir, 'mri_diff/mri_diff'))
      cmd = '%s %s %s' % (diffcmd, orig, ref)
      origname = rmext(op.basename(orig))
      cmd += ' --debug --diff diff-%s.mgz --log diff-%s.log' % (origname, origname)
      cmd += ' --thresh %f --res-thresh %f --geo-thresh %f' % (thresh, res_thresh, geo_thresh)
      if self.runcmd(cmd, fatal=False) != 0:
        errorExit('mri_diff of %s and %s failed' % (orig, ref))


  # overwrite a reference file with the observed output
  # this will re-tar and overwrite the primary testdata.tar.gz
  def regen(self, orig, ref):
    os.chdir(self.scriptdir)
    # create and enter a temporary dir to regenerate reference data
    tmp_dir = 'regeneration'
    rmdir(tmp_dir)
    os.makedirs(tmp_dir)
    os.chdir(tmp_dir)
    # extract the original testdata
    self.runcmd('tar -xzvf ../testdata.tar.gz')
    # replace the old with the new
    self.runcmd('mv -f %s %s' % (op.join(op.relpath(self.testdatadir), orig), ref))
    # tar up the testdata
    self.runcmd('tar -czvf testdata.tar.gz testdata/')
    os.chdir(self.scriptdir)
    # make sure the annex file is unlocked before replacing it
    self.runcmd('git annex unlock testdata.tar.gz')
    self.runcmd('mv -f %s/testdata.tar.gz testdata.tar.gz' % tmp_dir)
    rmdir(tmp_dir)


  # just delete the testdata directory and output any final messages
  def cleanup(self):
    self.cd(self.scriptdir)
    rmdir('testdata')
    # regeneration message
    if self.regenerate:
      print('testdata has been regenerated')
      print('make sure to run "git annex add testdata.tar.gz" to rehash before committing')
    else:
      print('%s[success]%s %s test has passed' % (term.green, term.end, op.basename(self.scriptdir)))
