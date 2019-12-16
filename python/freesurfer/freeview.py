import os
import tempfile
import numpy as np
import nibabel as nib

from . import error, run, Image, Overlay, Volume, Surface, collect_output


class Freeview:
    '''A visualization class that wraps a freeview command.

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
        self.volumes = []
        self.surfaces = []

    def vol(self, volume, **kwargs):
        '''Adds a new volume to the freeview command. If the volume provided is not a
        filepath, then the input will be saved as a volume in a temporary directory. 

        Args:
            volume: An existing volume filename, numpy array, or Volume object.
            name (str): Optional filename to apply to saved volumes. Do not include the extension.
            colormap (string): Colormap to apply to the volume.
            opacity (float): Opacity of the volume layer.
            opts (str): Additional options to append to the volume argument.
        '''

        # convert the input to a proper file (if it's not one already)
        filename = self._vol_to_file(volume)
        if filename is None:
            return
        self.volumes.append(filename + self._kwargs_to_flags(kwargs))

    def surf(self, surface, **kwargs):
        '''Adds a surface to the freeview session.

        Args:
            surface: An existing filename or Surface instance.
            name: Optional name for the loaded surface.
            overlay: A file, array, or Overlay instance to project onto the surface.
            mrisp: A file, array, or Image parameterization to project onto the surface.
            opts: Additional option string to append to the surface commandline argument.
        '''

        # retrieve the associated filename
        filename = self._surf_to_file(surface)
        if filename is None:
            return

        if kwargs.get('overlay'):
            kwargs['overlay'] = self._vol_to_file(kwargs['overlay'], force=Overlay)

        if kwargs.get('mrisp'):
            kwargs['mrisp'] = self._vol_to_file(kwargs['mrisp'] , force=Image)

        self.surfaces.append(filename + self._kwargs_to_flags(kwargs))

    def show(self, background=True, opts=''):
        '''Opens the configured freeview session.

        Args:
            background: Run freeview as a background process. Defaults to True.
            opts: Additional arguments to append to the command.
        '''

        # check for vgl on remote machinces since freeview can be a bit buggy
        cmd = self._vgl_wrapper() + 'freeview'

        # compile the rest of the command
        if self.volumes:
            cmd += ' -v ' + ' '.join(self.volumes)
        if self.surfaces:
            cmd += ' -f ' + ' '.join(self.surfaces)
        if opts:
            cmd += ' ' + opts

        # be sure to remove the temporary directory after freeview closes
        if self._tempdir is not None:
            cmd += ' ; rm -rf ' + self._tempdir

        # run it
        run(cmd, background=background)

    def _kwargs_to_flags(self, kwargs):

        if kwargs.get('name'):
            kwargs['name'] = kwargs['name'].replace(' ', '-')

        flags = kwargs.pop('opts', '')
        for key, value in kwargs.items():
            flags += ':%s=%s' % (key, str(value))

        return flags

    def _vol_to_file(self, volume, force=None):

        # if input is an already-existing file, no need to do anything
        if isinstance(volume, str):
            if not os.path.isfile(volume):
                error('volume file %s does not exist' % volume)
                return None
            return volume

        # if input is a numpy array, convert to the appropriate fs array type
        # and let the filename creation get handled below
        if isinstance(volume, np.ndarray):
            volume = self._convert_ndarray(volume) if force is None else force(volume.squeeze())
            if volume is None:
                error('cannot convert array of shape %s' % str(array.shape))
                return None

        # input is a Volume
        if isinstance(volume, Volume):
            filename = self._unique_filename('volume.mgz')
            volume.write(filename)
            return filename

        # input is an Image
        if isinstance(volume, Image):
            filename = self._unique_filename('image.mgz')
            volume.write(filename)
            return filename

        # input is an Overlay
        if isinstance(volume, Overlay):
            filename = self._unique_filename('overlay.mgz')
            volume.write(filename)
            return filename

        # as a final test, check if the input is possibly a nibabel image
        # we don't want nibabel to be required though, so ignore import errors
        try:
            import nibabel as nib
            if isinstance(volume, nib.spatialimages.SpatialImage):
                # convert from nib to fs volume, as it's easier to control the filename
                filename = self._unique_filename('volume.mgz')
                Volume(volume.get_data(), affine=volume.affine).write(filename)
                return filename
        except ImportError:
            pass

        error('invalid freeview volume type: %s' % type(volume))
        return None

    def _convert_ndarray(self, array):
        '''Converts a numpy array to the appropriate fs array instance.'''
        
        # remove any leading empty dimension
        if array.shape[0] == 1:
            array = array[0, ...]

        # remove and trailing empty dimension
        if array.shape[-1] == 1:
            array = array[..., 0]

        # compute number of base dims and frames
        # multi-frames are only assumed when array is 4D
        nframes = array.shape[-1] if array.ndim == 4 else 1
        array = array.squeeze()
        basedims = array.ndim if nframes == 1 else array.ndim - 1

        # construct corresponding array containers
        if basedims == 1:
            return Overlay(array)
        elif basedims == 2:
            return Image(array)
        elif basedims == 3:
            return Volume(array)
        else:
            return None

    def _surf_to_file(self, surface):
        
        # if input is an already-existing file, no need to do anything
        if isinstance(surface, str):
            if not os.path.isfile(surface):
                error('surface file %s does not exist' % surface)
                return None
            return surface
        
        # input is a Surface instance
        if isinstance(surface, Surface):
            filename = self._unique_filename('surface')
            surface.write(filename)
            return filename

        error('invalid freeview surface type: %s' % type(surface))
        return None

    def _get_temp_dir(self):
        # make the temporary directory if it does not already exist
        if self._tempdir is None:
            self._tempdir = tempfile.mkdtemp()
        return self._tempdir

    def _unique_filename(self, filename):
        '''Identifies a unique filename in the temporary directory.'''
        directory = self._get_temp_dir()
        # todoc
        fullpath = os.path.join(directory, filename)
        if not os.path.exists(fullpath):
            return fullpath
        # todoc
        name, ext = filename.split('.', 1)
        for n in range(2, 1000):
            fullpath = os.path.join(directory, '%s-%d.%s' % (name, n, ext))
            if not os.path.exists(fullpath):
                return fullpath

    def _vgl_wrapper(self):
        vgl = all([os.path.exists(path) for path in ('/etc/opt/VirtualGL/vgl_xauth_key', '/usr/pubsw/bin/vglrun')])
        local = any([os.environ.get('DISPLAY', '').endswith(string) for string in (':0', ':0.0')])
        if vgl and not local:
            if not 'NV-GLX' in collectOutput('xdpyinfo')[0]:
                return '/usr/pubsw/bin/vglrun '
        return ''


def fv(*args, **kwargs):
    '''Freeview wrapper to quickly load an arbitray number of volumes and surfaces. Inputs
    can be existing filenames, surfaces, volumes or numpy arrays. Use the `Freeview` class directly
    for a more specific configuration.

    Args:
        background: Run freeview as a background process. Defaults to True.
    '''
    background = kwargs.pop('background', True)
    opts = kwargs.pop('opts', '')

    fv = Freeview()
    for arg in args:

        # try to guess filetype if string
        if isinstance(arg, str):
            if arg.endswith(('.mgz', '.mgh', '.nii.gz', '.nii')):
                fv.vol(arg)
            elif arg.startswith(('lh.', 'rh.')) or arg.endswith('.stl'):
                fv.surf(arg)
            else:
                fv.vol(arg)

        # surface
        elif isinstance(arg, Surface):
            fv.surf(arg)

        # assume anything else is a volume
        else:
            fv.vol(arg)

    fv.show(background=background, opts=opts)


def fvoverlay(surface, overlay, background=True, opts=''):
    '''Freeview wrapper to quickly load an overlay onto a surface.

    Args:
        surface: An existing surface filename.
        overlay: An existing volume filename, a numpy array, or a nibabel image to apply as an overlay.
        background: Run freeview as a background process. Defaults to True.
    '''
    fv = Freeview()
    fv.surf(surface, overlay=overlay)
    fv.show(background=background, opts=opts)

