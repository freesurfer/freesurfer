import os
import tempfile
import numpy as np
import nibabel as nib

from . import error, run


class Freeview:
    '''A visualization class that wraps freeview.

    After configuring the wrapper, the freeview window can be opened with the show() method. Volumes and
    surfaces can be configured like so:

        fv = Freeview()
        fv.vol(volume)
        fv.vol(seg, colormap='lut', opacity=0.5)
        fv.surf(pialfile, overlay=thickness)
        fv.show()

    The name, colormap, and opacity options for vol() are just easy ways to specify the key-value argument
    options on the commandline. Additional options can be applied to a volume or surface with the 'opts' option.
    So, for example, the following:

        fv = Freeview()
        fv.vol(volume, name='seg1', colormap='lut', opacity=0.5)
        fv.show()

    Is the same as:

        fv = Freeview()
        fv.vol(volume, opts=':name=seg1:colormap=lut:opacity=0.5')
        fv.show()

    Additional flags can be passed to freeview in the show() method:

        fv = Freeview()
        fv.show(opts='--slice 5 5 5 --verbose')

    For a quicker but more limited way to view volumes or overlays, see `fv()` and `fvoverlay()`.
    '''
    def __init__(self):
        self._tempdir = None
        self.items = {'volumes': [], 'surfaces': []}

    def vol(self, volume, affine=np.eye(4), name=None, colormap=None, opacity=None, opts=''):
        '''Adds a new volume to the freeview command. If the volume provided is not a
        filepath, then the input will be saved as a volume in a temporary directory. 

        Args:
            volume: An existing volume filename, a numpy array, or a nibabel image.
            name (str): Optional filename to apply to saved volumes. Do not include the extension.
            affine (numpy.ndarray): A 4x4 affine transform to apply to the saved volume. The identity matrix 
                is used by default.
            colormap (string): Colormap to apply to the volume.
            opacity (float): Opacity of the volume layer.
            opts (str): Additional options to append to the volume argument.
        '''
        filename = self._makeVolume(volume, 'volume', affine)
        if filename is None:
            return

        if name is not None:
            opts += ':name=' + name
        if colormap is not None:
            opts += ':colormap=' + colormap
        if opacity is not None:
            opts += ':opacity=' + str(opacity)

        self.items['volumes'].append(filename + opts)

    def surf(self, surface, overlay=None, mrisp=None, opts=''):
        '''Adds a new surface to the freeview command. Currently, only existing surface file
        must be provided as input.

        Args:
            surface: An existing surface filename.
            overlay: An existing volume filename, a numpy array, or a nibabel image to apply as an overlay.
            mrisp: An existing volume filename, a numpy array, or a nibabel image to apply as an mrisp.
            opts (str): Additional options to append to the surface argument.
        '''
        if isinstance(surface, str):
            # input is an already-existing surface file
            if not os.path.exists(surface):
                error('surface file %s does not exist' % surface)
                return
            filename = surface
        else:
            error('invalid freeview surface argument type %s' % type(surface))
            return

        if overlay is not None:
            opts += ':overlay=' + self._makeVolume(overlay, 'overlay', np.eye(4))
        if mrisp is not None:
            opts += ':mrisp=' + self._makeVolume(mrisp, 'mrisp', np.eye(4))

        self.items['surfaces'].append(filename + opts)

    def show(self, background=True, opts=''):
        '''Runs the configured freeview command.

        Args:
            background: Run freeview as a background process. Defaults to True.
            opts: Additional arguments to append to the command.
        '''

        cmd = 'freeview'

        # use vgl if remote since freeview can be a bit buggy
        vgl = all([os.path.exists(path) for path in ('/etc/opt/VirtualGL/vgl_xauth_key', '/usr/pubsw/bin/vglrun')])
        local = any([os.environ.get('DISPLAY', '').endswith(string) for string in (':0', ':0.0')])
        if vgl and not local:
            cmd = '/usr/pubsw/bin/vglrun ' + cmd

        if self.items['volumes']:
            cmd += ' -v ' + ' '.join(self.items['volumes'])
        if self.items['surfaces']:
            cmd += ' -f ' + ' '.join(self.items['surfaces'])
        if opts:
            cmd += ' ' + opts
        if self._tempdir is not None:
            # delete the temporary directory after freeview closes
            cmd += ' ; rm -rf ' + self._tempdir

        run(cmd, background=background)

    def _makeVolume(self, volume, name, affine):
        if isinstance(volume, str):
            # input is an already-existing volume file
            if not os.path.exists(volume):
                error('volume file %s does not exist' % volume)
                return None
            filename = volume
        elif isinstance(volume, nib.spatialimages.SpatialImage):
            # input is a nibabel image
            volume.set_filename(os.path.join(self._makeTempDir(), name))
            volume.set_filename(self._validName(volume.get_filename()))
            nib.save(volume, volume.get_filename())
            filename = volume.get_filename()
        elif isinstance(volume, np.ndarray):
            # input is a nifty array
            filename = self._validName(os.path.join(self._makeTempDir(), name + '.nii.gz'))
            nib.save(nib.Nifti1Image(volume, affine), filename)
        else:
            error('invalid freeview volume argument type %s' % type(volume))
            return None
        return filename

    def _makeTempDir(self):
        # make the temporary directory if it does not already exist
        if self._tempdir is None:
            self._tempdir = tempfile.mkdtemp()
        return self._tempdir

    def _validName(self, filename):
        # identifies a unique filename in the temporary directory
        if not os.path.exists(filename):
            return filename
        else:
            directory, file = os.path.split(filename)
            for n in range(2, 1000):
                name, ext = file.split('.', 1)
                unique = os.path.join(directory, '%s%d.%s' % (name, n, ext))
                if not os.path.exists(unique):
                    break
            return unique


def fv(*args, **kwargs):
    '''Freeview wrapper to quickly load an arbitray number of volumes. Inputs can be
    existing volume filenames, numpy arrays, or nibabel images.

    Args:
        background: Run freeview as a background process. Defaults to True.
    '''
    background = kwargs.pop('background', True)
    opts = kwargs.pop('opts', '')

    freeview = Freeview()
    for arg in args:
        freeview.vol(arg)
    freeview.show(background=background, opts=opts)


def fvoverlay(surface, overlay, background=True, opts=''):
    '''Freeview wrapper to quickly load an overlay onto a surface.

    Args:
        surface: An existing surface filename.
        overlay: An existing volume filename, a numpy array, or a nibabel image to apply as an overlay.
        background: Run freeview as a background process. Defaults to True.
    '''
    freeview = Freeview()
    freeview.surf(surface, overlay=overlay)
    freeview.show(background=background, opts=opts)
