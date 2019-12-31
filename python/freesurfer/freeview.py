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

    For a quicker but more limited way to view volumes or overlays, see `fv()` and `fvoverlay()`.
    '''
    def __init__(self):
        self.tempdir = None
        self.flags = []

    def vol(self, volume, **kwargs):
        '''
        Loads a volume in the sessions. If the volume provided is not a filepath,
        then the input will be saved as a volume in a temporary directory. Any
        key/value tags allowed on the command line can be provided as arguments.

        Args:
            volume: An existing volume filename, numpy array, or fs array container.
            opts: Additional options to append to the volume argument.
        '''

        # convert the input to a proper file (if it's not one already)
        filename = self._vol_to_file(volume)
        if filename is None:
            return

        # build the volume flag
        flag = '-v ' + filename + self._kwargs_to_tags(kwargs)
        self.add_flag(flag)

    def surf(self, surface, overlay=None, mrisp=None, sphere=None, **kwargs):
        '''
        Loads a surface in the freeview session. If the surface provided is not
        a filepath, then the input will be saved in a temporary directory. Any
        key/value tags allowed on the command line can be provided as arguments.

        Args:
            surface: An existing filename or Surface instance.
            overlay: A file, array, or Overlay instance to project onto the surface.
            mrisp: A file, array, or Image parameterization to project onto the surface.
            opts: Additional option string to append to the surface commandline argument.
        '''

        # convert the input to a proper file (if it's not one already)
        filename = self._surf_to_file(surface)
        if filename is None:
            return

        # if overlay is provided as an array, make sure it's converted
        if overlay is not None:
            kwargs['overlay'] = self._vol_to_file(overlay, force=Overlay)

        # if mrisp is provided as an array, make sure it's converted
        if mrisp is not None:
            kwargs['mrisp'] = self._vol_to_file(mrisp, force=Image)

        # if sphere is provided as an array, make sure it's converted
        if sphere is not None:
            kwargs['sphere'] = self._surf_to_file(sphere)

        # build the surface flag
        flag = '-f ' + filename + self._kwargs_to_tags(kwargs)
        self.add_flag(flag)

    def show(self, background=True, opts=''):
        '''Opens the configured freeview session.

        Args:
            background: Run freeview as a background process. Defaults to True.
            opts: Additional arguments to append to the command.
        '''

        # compile the command
        command = '%s freeview %s %s' % (self._vgl_wrapper(), opts, ' '.join(self.flags))

        # be sure to remove the temporary directory (if it exists) after freeview closes
        if self.tempdir:
            command = '%s ; rm -rf %s' % (command, self.tempdir)

        # run it
        run(command, background=background)

    def add_flag(self, flag):
        '''Adds a flag to the freeview command.'''
        self.flags.append(flag)

    def _kwargs_to_tags(self, kwargs):
        '''Converts a kwargs dictionary to a freeview key/value tag string.'''

        # make sure there are no spaces in the "name" field
        if kwargs.get('name'):
            kwargs['name'] = kwargs['name'].replace(' ', '-')

        # opts is reserved for hardcoded tags
        tags = kwargs.pop('opts', '')
        for key, value in kwargs.items():
            tags += ':%s=%s' % (key, str(value))

        return tags

    def _vol_to_file(self, volume, force=None):
        '''
        Converts an unknown volume type (whether it's a filename, array, or
        other object) into a valid file.
        '''

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
        '''
        Converts a numpy array to the appropriate fs array instance.
        The appropriate type is guessed based on the shape of the array.
        '''
        
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
        '''
        Converts an unknown surface type (whether it's a filename or
        actual surface) into a valid file.
        '''
        
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
        '''Returns the temporary directory. If it does not exist, it gets created.'''
        if self.tempdir is None:
            self.tempdir = tempfile.mkdtemp()
        return self.tempdir

    def _unique_filename(self, filename):
        '''Identifies a unique filename in the temporary directory.'''

        directory = self._get_temp_dir()
        
        # check if it's unique
        fullpath = os.path.join(directory, filename)
        if not os.path.exists(fullpath):
            return fullpath

        # append numbers until a unique filename is created (stop after 10,000 tries)
        name_ext = filename.split('.', 1)
        name = name_ext[0]
        ext = '.' + name_ext[1] if len(name_ext) == 2 else ''
        for n in range(2, 10000):
            fullpath = os.path.join(directory, '%s-%d%s' % (name, n, ext))
            if not os.path.exists(fullpath):
                return fullpath
        raise RuntimeError('could not generate a unique filename for "%s" after trying many times' % (filename))

    def _vgl_wrapper(self):
        '''
        Freeview can be buggy when run remotely. This is a martinos-specific solution
        that wraps the process with VGL.
        '''
        vgl = all([os.path.exists(path) for path in ('/etc/opt/VirtualGL/vgl_xauth_key', '/usr/pubsw/bin/vglrun')])
        local = any([os.environ.get('DISPLAY', '').endswith(string) for string in (':0', ':0.0')])
        if vgl and not local:
            if not 'NV-GLX' in collect_output('xdpyinfo')[0]:
                return '/usr/pubsw/bin/vglrun '
        return ''


def fv(*args, **kwargs):
    '''
    Freeview wrapper to quickly load an arbitray number of volumes and surfaces. Inputs
    can be existing filenames, surfaces, volumes or numpy arrays. Use the `Freeview` class directly
    to configure a more detailed session.

    Args:
        opts: Additional string of flags to add to the command.
        background: Run freeview as a background process. Defaults to True.
    '''
    background = kwargs.pop('background', True)
    opts = kwargs.pop('opts', '')

    fv = Freeview()
    for arg in args:
        if isinstance(arg, str):
            # try to guess filetype if string
            if arg.endswith(('.mgz', '.mgh', '.nii.gz', '.nii')):
                fv.vol(arg)
            elif arg.startswith(('lh.', 'rh.')) or arg.endswith('.stl'):
                fv.surf(arg)
            else:
                fv.vol(arg)
        elif isinstance(arg, Surface):
            # surface
            fv.surf(arg)
        else:
            # assume anything else is a volume
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
