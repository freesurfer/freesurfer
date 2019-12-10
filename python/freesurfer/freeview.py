import os
import tempfile
import numpy as np
import nibabel as nib

from . import error, run, Volume, Surface, collectOutput


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

    def vol(self, volume, affine=None, name=None, colormap=None, opacity=None, opts=''):
        '''Adds a new volume to the freeview command. If the volume provided is not a
        filepath, then the input will be saved as a volume in a temporary directory. 

        Args:
            volume: An existing volume filename, numpy array, or Volume object.
            name (str): Optional filename to apply to saved volumes. Do not include the extension.
            affine (numpy.ndarray): A 4x4 affine transform to apply to the saved volume.
            colormap (string): Colormap to apply to the volume.
            opacity (float): Opacity of the volume layer.
            opts (str): Additional options to append to the volume argument.
        '''
        filename = self._get_volume_file(volume, 'volume', affine)
        if filename is None:
            return

        if name is not None:
            opts += ':name=' + name.replace(' ', '-')
        if colormap is not None:
            opts += ':colormap=' + colormap
        if opacity is not None:
            opts += ':opacity=' + str(opacity)

        self.items['volumes'].append(filename + opts)

    def surf(self, surface, name=None, overlay=None, mrisp=None, opts=''):
        '''Adds a new surface to the freeview command.

        Args:
            surface: An existing surface filename or Surface object.
            name (str): Optional filename to apply to saved surfaces.
            overlay: An existing volume filename, a numpy array, or a nibabel image to apply as an overlay.
            mrisp: An existing volume filename, a numpy array, or a nibabel image to apply as an mrisp.
            opts (str): Additional options to append to the surface argument.
        '''
        filename = self._get_surface_file(surface, 'surface')
        if filename is None:
            return

        if name is not None:
            opts += ':name=' + name.replace(' ', '-')
        if overlay is not None:
            overlay = self._overlay_to_vol(overlay)
            opts += ':overlay=' + self._get_volume_file(overlay, 'overlay', np.eye(4))
        if mrisp is not None:
            opts += ':mrisp=' + self._get_volume_file(mrisp, 'mrisp', np.eye(4))

        self.items['surfaces'].append(filename + opts)

    def show(self, background=True, opts=''):
        '''Runs the configured freeview command.

        Args:
            background: Run freeview as a background process. Defaults to True.
            opts: Additional arguments to append to the command.
        '''

        if not self.items['volumes'] and not self.items['surfaces']:
            error('nothing to load in freeview - not opening')
            return

        cmd = 'freeview'

        # use vgl if remote since freeview can be a bit buggy
        vgl = all([os.path.exists(path) for path in ('/etc/opt/VirtualGL/vgl_xauth_key', '/usr/pubsw/bin/vglrun')])
        local = any([os.environ.get('DISPLAY', '').endswith(string) for string in (':0', ':0.0')])
        if vgl and not local:
            if not 'NV-GLX' in collectOutput('xdpyinfo')[0]:
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

    def _get_volume_file(self, volume, name, affine):
        if isinstance(volume, str):
            # input is an already-existing volume file
            if not os.path.exists(volume):
                error('volume file %s does not exist' % volume)
                return None
            filename = volume
        elif isinstance(volume, nib.spatialimages.SpatialImage):
            # input is a nibabel image
            volume.set_filename(os.path.join(self._get_temp_dir(), name))
            volume.set_filename(self._get_valid_name(volume.get_filename()))
            nib.save(volume, volume.get_filename())
            filename = volume.get_filename()
        elif isinstance(volume, Volume):
            filename = self._get_valid_name(os.path.join(self._get_temp_dir(), name + '.mgz'))
            volume.write(filename)
        elif isinstance(volume, np.ndarray):
            # int64 MRI IO isn't very stable, so convert to int32 for now
            if volume.dtype == 'int64':
                volume = volume.astype('int32')
            # make sure dimensions are valid
            if volume.ndim > 4:
                volume = volume.squeeze()
                if volume.ndim > 4:
                    error('freeview input array has %d dimensions' % volume.ndim)
                    return None
            if volume.ndim < 3:
                volume = volume[..., np.newaxis]
            # input is a nifty array
            filename = self._get_valid_name(os.path.join(self._get_temp_dir(), name + '.mgz'))
            vol = Volume(volume)
            if affine is not None:
                vol.affine = affine
            vol.write(filename)
        else:
            error('invalid freeview volume argument type %s' % type(volume))
            return None
        return filename

    def _get_surface_file(self, surface, name):
        if isinstance(surface, str):
            # input is an already-existing surface file
            if not os.path.exists(surface):
                error('surface file %s does not exist' % surface)
                return None
            filename = surface
        elif isinstance(surface, Surface):
            filename = self._get_valid_name(os.path.join(self._get_temp_dir(), name))
            surface.write(filename)
        else:
            error('invalid freeview surface argument type %s' % type(surface))
            return None
        return filename

    def _get_temp_dir(self):
        # make the temporary directory if it does not already exist
        if self._tempdir is None:
            self._tempdir = tempfile.mkdtemp()
        return self._tempdir

    def _get_valid_name(self, filename):
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

    def _overlay_to_vol(self, overlay):
        if not isinstance(overlay, np.ndarray):
            return overlay
        if len(overlay.shape) == 1:
            return overlay.reshape((overlay.shape[0], 1, 1, 1))
        elif len(overlay.shape) == 2:
            return overlay.reshape((overlay.shape[0], 1, 1, overlay.shape[1]))
        else:
            error('overlay cannot be 3D!') 



def fv(*args, **kwargs):
    '''Freeview wrapper to quickly load an arbitray number of volumes and surfaces. Inputs
    can be existing filenames, surfaces, volumes or numpy arrays. Use the `Freeview` class directly
    for a more specific configuration.

    Args:
        background: Run freeview as a background process. Defaults to True.
    '''
    background = kwargs.pop('background', True)
    opts = kwargs.pop('opts', '')

    freeview = Freeview()
    for arg in args:
        if isinstance(arg, Surface):
            freeview.surf(arg)
        else:
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

