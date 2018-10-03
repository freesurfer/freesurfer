import os
import tempfile
import numpy as np
import nibabel as nib

from . import error, run


def fv(*args, **kwargs):
    """Freeview wrapper that accepts filenames, nibabel images, and numpy arrays.

    Optional Args:
        affine: affine transform for saving numpy arrays (default: identity matrix).
        background: run freeview as a background process (default: True).

    """
    affine = kwargs.pop('affine', np.eye(4))
    background = kwargs.pop('background', True)

    surfaces = []
    volumes = []
    saved_volumes = 0
    tmpdir = ''

    for arg in args:

        if isinstance(arg, str):
            # -- string --
            if not os.path.exists(arg):
                error("file '%s' does not exist" % arg)
                continue
            filename = os.path.basename(arg)
            # volume
            if filename.endswith(('.mgh', '.mgz', '.nii', '.nii.gz')):
                volumes.append(arg)
            # surface
            elif filename.startswith(('lh.', 'rh.')):
                surfaces.append(arg)
            # unknown
            else:
                error("could not determine file type of '%s'" % arg)

        elif isinstance(arg, nib.spatialimages.SpatialImage):
            # -- nibabel image -
            tmpdir = makeTmpDir(tmpdir)
            saved_volumes += 1
            arg.set_filename(os.path.join(tmpdir, 'vol%d' % saved_volumes))
            nib.save(arg, arg.get_filename())
            volumes.append(arg.get_filename())

        elif isinstance(arg, np.ndarray):
            # -- numpy array --
            tmpdir = makeTmpDir(tmpdir)
            saved_volumes += 1
            fname = os.path.join(tmpdir, 'vol%d.nii' % saved_volumes)
            nib.save(nib.Nifti1Image(arg, affine), fname)
            volumes.append(fname)

        else:
            # -- invalid --
            error('invalid fv argument type %s' % type(arg))

    # freeview command
    cmd = 'freeview'
    if volumes:  cmd += ' -v ' + ' '.join(volumes)
    if surfaces: cmd += ' -f ' + ' '.join(surfaces)
    # make sure tmp data gets cleaned up after freeview closes
    if tmpdir: cmd += ' && rm -rf %s' % tmpdir
    run(cmd, background=background)


def makeTmpDir(tmpdir):
    if not tmpdir:
        tmpdir = tempfile.mkdtemp()
    return tmpdir
