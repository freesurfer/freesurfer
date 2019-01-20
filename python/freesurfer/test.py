import os, os.path as op
import sys
import shutil

from . import run, term, errorExit, ArgParser


class RegressionTest:
    '''
    This class is built specifically for configuring regression tests for
    freesurfer binaries in the main source tree. This assumes all test data is
    stored as a 'testdata.tar.gz' tarball in the directory where the test script
    lives. Commands should be tested using run(), and outputs should be compared
    against reference data by using the diff(), mridiff, or surfdiff() functions.
    Example:

      rt = fst.RegressionTest()

      rt.run('mri_normalize -mprage nu.mgz T1.mgz')
      rt.mridiff('T1.mgz', 'T1.ref.mgz')

      rt.run('mri_normalize -gentle orig.mgz gentle.mgz')
      rt.mridiff('gentle.mgz', 'gentle.ref.mgz')

      rt.cleanup()

    The initialization of a RegressionTest object will parse the command line and
    look for the --regenerate flag, which will overwrite the reference data with any
    produced outputs specified in the diff functions.
    '''

    def __init__(self):
        # parse command line arguments
        parser = ArgParser()
        parser.add_argument('--regenerate', action='store_true', help='Regenerate the reference data.')
        parser.add_argument('--keep-data', action='store_true', help='Keep the testdata dir even after success.')
        args = parser.parse_args()

        # get directory of the test script
        self.keepdata = args.keep_data
        self.testdir = os.getcwd()
        self.testdatadir = op.join(self.testdir, 'testdata')
        self.scriptdir = op.dirname(op.realpath(sys.argv[0]))
        self.testdatatar = op.join(self.scriptdir, 'testdata.tar.gz')
        self.fshome = self.findPath(self.scriptdir, 'distribution')
        self._setEnvironment()

        # if regenerating testdata...
        self.regenerate = args.regenerate
        if self.regenerate:
            print('regenerating %s testdata' % op.basename(self.scriptdir))
            # make the tmp generation directory
            self.regeneration_dir = op.join(self.testdir, 'testdata_regeneration')
            shutil.rmtree(self.regeneration_dir, ignore_errors=True)
            os.makedirs(self.regeneration_dir)
            os.chdir(self.regeneration_dir)
            # extract the original testdata to get replaced later
            self._runcmd('tar -xzvf "%s"' % self.testdatatar)
            os.chdir(self.testdir)

    def _setEnvironment(self):
        # setup the default testing environment
        os.environ['FREESURFER_HOME'] = self.fshome
        os.environ['SUBJECTS_DIR'] = self.testdatadir
        os.environ['FSLOUTPUTTYPE'] = 'NIFTI_GZ'
        # setup martinos license if needed
        haveLicense = os.environ.get('FS_LICENSE') or \
                      op.exists(op.join(self.fshome, '.license')) or \
                      op.exists(op.join(self.fshome, 'license.txt'))
        if not haveLicense and op.exists('/space/freesurfer/.license'):
            os.environ["FS_LICENSE"] = '/space/freesurfer/.license'

    # recursively search upwards for a target
    # this is used to locate the 'distribution' folder or the compiled 'mri_diff'
    # in the build directory
    def findPath(self, startpath, target):
        distpath = op.join(startpath, target)
        if not op.exists(distpath):
            parent = op.abspath(op.join(startpath, '..'))
            if parent == '/':
                errorExit("max recursion reached - could not locate '%s'" % target)
            return self.findPath(parent, target)
        return distpath

    # print and run a command, and exit if it fails
    def _runcmd(self, cmd, fatal=True):
        print(term.bold('>> ' + cmd))
        ret = run(cmd)
        if fatal and ret != 0:
            errorExit('"%s" failed' % cmd)
        return ret

    def run(self, cmd, threads=None, allow_failure=False, expect_failure=False):
        '''Runs a command (this is what should be called from the actual test script).'''
        os.chdir(self.testdir)
        # extract the testdata
        shutil.rmtree(self.testdatadir, ignore_errors=True)
        self._runcmd('tar -xzvf "%s"' % self.testdatatar)
        os.chdir(self.testdatadir)
        # set number of OMP threads if necessary
        if threads is not None:
            os.environ['OMP_NUM_THREADS'] = str(threads)
            print('setting OMP_NUM_THREADS = %s' % os.environ['OMP_NUM_THREADS'])
        # run the test command - assumes the command is located in the testdir (above testdatadir)
        ret = self._runcmd('../' + cmd, fatal=False)
        if expect_failure:
            if ret == 0:
                errorExit('test command "%s" returned 0, but expected a failure')
        elif ret != 0 and not allow_failure:
            errorExit('test command "%s" failed' % cmd)

    def diff(self, orig, ref, ignore_comments=False):
        '''Runs the standard unix diff command.'''
        if self.regenerate:
            self._regen(orig, ref)
        else:
            os.chdir(self.testdatadir)
            cmd = 'diff %s %s' % (orig, ref)
            if ignore_comments:
                cmd += " -I '#'"
            if self._runcmd(cmd, fatal=False) != 0:
                errorExit('diff of %s and %s failed' % (orig, ref))

    def mridiff(self, orig, ref, thresh=0.0, res_thresh=1e-6, geo_thresh=8e-6, flags=""):
        '''Runs a diff on two volumes (calls mri_diff, which must be already built in the source directory).'''
        if self.regenerate:
            self._regen(orig, ref)
        else:
            os.chdir(self.testdatadir)
            diffcmd = op.relpath(self.findPath(self.testdatadir, 'mri_diff/mri_diff'))
            cmd = '%s %s %s' % (diffcmd, orig, ref)
            origname = op.basename(orig)
            for ext in ('.mgh', '.mgz', '.nii', '.nii.gz'):
                if orig.endswith(ext):
                    origname = origname[:-len(ext)]
            cmd += ' --debug --diff diff-%s.mgz --log diff-%s.log' % (origname, origname)
            cmd += ' --thresh %f --res-thresh %f --geo-thresh %f %s' % (thresh, res_thresh, geo_thresh, flags)
            if self._runcmd(cmd, fatal=False) != 0:
                errorExit('mri_diff of %s and %s failed' % (orig, ref))

    def surfdiff(self, orig, ref, flags=""):
        '''Runs a diff on two surfs (calls `mris_diff`, which must be already built in the source directory).'''
        if self.regenerate:
            self._regen(orig, ref)
        else:
            os.chdir(self.testdatadir)
            diffcmd = op.relpath(self.findPath(self.testdatadir, 'mris_diff/mris_diff'))
            cmd = '%s %s %s --debug %s' % (diffcmd, orig, ref, flags)
            if self._runcmd(cmd, fatal=False) != 0:
                errorExit('mris_diff of %s and %s failed' % (orig, ref))
            os.chdir(self.testdir)  # return to testdir

    def _regen(self, orig, ref):
        # overwrites a reference file with the observed output
        # this will re-tar and overwrite the primary testdata.tar.gz
        os.chdir(self.regeneration_dir)
        # replace the old with the new
        self._runcmd('mv -f %s %s' % (op.join(op.relpath(self.testdatadir), orig), op.join('testdata', ref)))
        os.chdir(self.testdatadir)

    def cleanup(self):
        '''Deletes the testdata directory and outputs any final messages. This should only get called when the
        whole test finishes successfully.'''
        os.chdir(self.testdir)
        if not self.keepdata:
            shutil.rmtree('testdata', ignore_errors=True)
        if self.regenerate:
            # tar the regenerated data
            os.chdir(self.regeneration_dir)
            self._runcmd('tar -czvf testdata.tar.gz testdata')
            # make sure the annex file is unlocked before replacing it
            os.chdir(self.scriptdir)
            self._runcmd('git annex unlock testdata.tar.gz')
            self._runcmd('mv -f %s .' % op.join(self.regeneration_dir, 'testdata.tar.gz'))
            shutil.rmtree(self.regeneration_dir)
            print('testdata has been regenerated')
            print('make sure to run "git annex add testdata.tar.gz" to rehash before committing')
            os.chdir(self.testdir)
        else:
            print('%s %s test has passed' % (term.green('[success]'), op.basename(self.scriptdir)))
