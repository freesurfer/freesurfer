import os
import tempfile
import numpy as np
import nibabel as nib

from . import error, run


def fv(*args, **kwargs):
    """Freeview wrapper that accepts filenames, nibabel images, and numpy arrays.

    Optional Args:
        affine: 4x4 affine transform for saving numpy arrays. Uses identity matrix by default.
        mrisp: Parameterization to load on top of the first surface file provided.
        overlay: Overlay to load on top of the first surface file provided.
        background (bool): Run freeview as a background process. Defaults to True.
    """
    affine = kwargs.pop('affine', np.eye(4))
    mrisp = kwargs.pop('mrisp', None)
    overlay = kwargs.pop('overlay', None)
    background = kwargs.pop('background', True)

    # create directory to save temporary volumes and surfaces
    tmpdir = tempfile.mkdtemp()

    surfaces = []
    volumes = []
    for arg in args:
        if isinstance(arg, str):
            # assume all string inputs are existing filenames
            if not os.path.exists(arg):
                error('file %s does not exist' % arg)
                continue
            filename = os.path.basename(arg)
            if filename.endswith(('.mgh', '.mgz', '.nii', '.nii.gz')):
                volumes.append(arg)
            elif filename.startswith(('lh.', 'rh.')):
                surfaces.append(arg)
            else:
                error('cannot determine filetype of %s' % arg)
        else:
            # otherwise, assume input is a volume
            volumes.append(_mkvolume(arg, os.path.join(tmpdir, 'vol%d' % (len(volumes)+1)), affine=affine))

    if mrisp is not None:
        if not surfaces:
            error('cannot load mrisp if no surface is provided')
        surfaces[0] += ':mrisp=%s' % _mkvolume(mrisp, os.path.join(tmpdir, 'mrisp'))

    if overlay is not None:
        if not surfaces:
            error('cannot load overlay if no surface is provided')
        surfaces[0] += ':overlay=%s' % _mkvolume(overlay, os.path.join(tmpdir, 'overlay'))

    # formulate the command
    cmd = 'freeview'
    # use martinos vgl if running remotely
    display = os.environ.get('DISPLAY', '')
    if os.path.exists('/etc/opt/VirtualGL/vgl_xauth_key') and \
        not display.endswith(':0') and not display.endswith(':0.0'):
        cmd = '/usr/pubsw/bin/vglrun freeview'

    if volumes:
        cmd += ' -v ' + ' '.join(volumes)
    if surfaces:
        cmd += ' -f ' + ' '.join(surfaces)

    cmd += ' ; rm -rf ' + tmpdir
    run(cmd, background=background)


def _mkvolume(arg, filename, affine=np.eye(4)):
    if isinstance(arg, str):
        return arg
    elif isinstance(arg, nib.spatialimages.SpatialImage):
        arg.set_filename(filename)
        nib.save(arg, arg.get_filename())
        return arg.get_filename()
    elif isinstance(arg, np.ndarray):
        path = '%s.nii' % filename
        nib.save(nib.Nifti1Image(arg, affine), path)
        return path
    else:
        error('invalid fv argument type %s' % type(arg))
