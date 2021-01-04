import os
import tempfile
import numpy as np

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

    class OverlayTag:
        '''Configuration for overlay tags. See surf() for usage.'''
        def __init__(self, data, name=None, threshold=None, opacity=None):
            self.data = data
            self.name = name
            self.threshold = threshold
            self.opacity = opacity

    class AnnotTag:
        '''Configuration for annotation tags. See surf() for usage.'''
        def __init__(self, data, lut=None, name=None):
            self.data = Overlay(data.squeeze()) if isinstance(data, np.ndarray) else data
            if lut is not None:
                self.data.lut = lut
            self.name = name

    class MRISPTag:
        '''Configuration for mrisp tags. See surf() for usage.'''
        def __init__(self, data, name=None):
            self.data = data
            self.name = name

    def __init__(self):
        self.tempdir = None
        self.flags = []

    def copy(self):
        '''
        Returns a copy of the current freeview configuration. This allows for multiple
        freeview windows loading the same base file (without resaving it).
        '''
        copied = Freeview()
        copied.flags = self.flags.copy()
        copied.tempdir = self.tempdir
        return copied

    def vol(self, volume, swap_batch_dim=False, **kwargs):
        '''
        Loads a volume in the sessions. If the volume provided is not a filepath,
        then the input will be saved as a volume in a temporary directory. Any
        key/value tags allowed on the command line can be provided as arguments.

        Args:
            volume: An existing volume filename, numpy array, or fs array container.
            swap_batch_dim: Move the first axis to the last if input is a numpy array. Default is False.
            opts: Additional options to append to the volume argument.
        '''

        # convert the input to a proper file (if it's not one already)
        filename = self._vol_to_file(volume, swap_batch_dim=swap_batch_dim)
        if filename is None:
            return

        # build the volume flag
        flag = '-v ' + filename + self._kwargs_to_tags(kwargs)
        self.add_flag(flag)

    def surf(self, surface, overlay=None, annot=None, mrisp=None, sphere=None, curvature=None, **kwargs):
        '''
        Loads a surface in the freeview session. If the surface provided is not
        a filepath, then the input will be saved in a temporary directory. Any
        key/value tags allowed on the command line can be provided as arguments.

        Overlays can be provided as filenames, numpy arrays, or Overlay instances, but in order
        to configure overlays with custom visualization options, use the OverlayTag class to specify
        things like desired filename and threshold:

            overlay = fs.Freeview.OverlayTag(thickness, name='thickness', threshold=(1, 3))
            fv.surf(surface, overlay=overlay)

        The first argument to the OverlayTag constructor can be a filename, numpy array, or Overlay
        instance. A list of multiple OverlayTags can be provided as input to the overlay parameter as
        well. A similar configuration class exists for mrisps, named MRISPTag.

        Args:
            surface: An existing filename or Surface instance.
            overlay: A file, array, Overlay, or OverlayTag instance to project onto the surface. Multiple overlays can
                be provided with a list.
            annot: An Overlay or AnnotTag instance to annotate the surface. Multiple
                annotations can be provided with a list. Overlays must have embedded lookup tables.
            mrisp: A file, array, Image, or MRISPTag to project onto the surface. Multiple parameterizations can
                be provided with a list.
            curvature: A file, array, or Overlay instance to load as the surface curvature.
            opts: Additional option string to append to the surface commandline argument.
        '''

        # convert the input to a proper file (if it's not one already)
        filename = self._surf_to_file(surface)
        if filename is None:
            return

        # if curvature is provided as a np array, make sure it's converted
        if curvature is not None:
            kwargs['curvature'] = self._vol_to_file(curvature, force=Overlay)

        # configure (potentially multiple) overlays
        if overlay is not None:
            overlay = list(overlay) if isinstance(overlay, (list, tuple)) else [overlay]
            for ol in overlay:
                config = ol if isinstance(ol, Freeview.OverlayTag) else Freeview.OverlayTag(ol)
                tag = ':overlay=%s' % self._vol_to_file(config.data, name=config.name, force=Overlay)
                if config.threshold is not None:
                    tag += ':overlay_threshold=%s' % (','.join(str(x) for x in config.threshold))
                if config.opacity is not None:
                    tag += ':overlay_opacity=%f' % config.opacity
                kwargs['opts'] = tag + kwargs.get('opts', '')

        # configure (potentially multiple) annots
        if annot is not None:
            annot = list(annot) if isinstance(annot, (list, tuple)) else [annot]
            for an in annot:
                config = an if isinstance(an, Freeview.AnnotTag) else Freeview.AnnotTag(an)
                if config.name is None:
                    config.name = 'annotation'
                annot_filename = self._vol_to_file(config.data, name=config.name, force=Overlay, ext='annot')
                if annot_filename is not None:
                    tag = ':annot=%s' % annot_filename
                    kwargs['opts'] = tag + kwargs.get('opts', '')

        # configure (potentially multiple) mrisps
        if mrisp is not None:
            mrisp = list(mrisp) if isinstance(mrisp, (list, tuple)) else [mrisp]
            for sp in mrisp:
                config = sp if isinstance(sp, Freeview.MRISPTag) else Freeview.MRISPTag(sp)
                tag = ':mrisp=%s' % self._vol_to_file(config.data, name=config.name, force=Image)
                kwargs['opts'] = tag + kwargs.get('opts', '')

        # if sphere is provided as an array, make sure it's converted
        if sphere is not None:
            kwargs['sphere'] = self._surf_to_file(sphere)

        # build the surface flag
        flag = '-f ' + filename + self._kwargs_to_tags(kwargs)
        self.add_flag(flag)

    def show(self, background=True, title=None, opts='', verbose=False, noclean=False, threads=None):
        '''Opens the configured freeview session.

        Args:
            background: Run freeview as a background process. Defaults to True.
            title: Title for the freeview window. Defaults to None.
            opts: Additional arguments to append to the command.
            verbose: Print the freeview command before running. Defaults to False.
            noclean: Do not remove temporary directory for debugging purposes.
            threads: Set number of OMP threads available to freeview.
        '''

        # compile the command
        command = '%s freeview %s %s' % (self._vgl_wrapper(), opts, ' '.join(self.flags))

        if title is not None:
            command += ' -subtitle "%s"' % title.replace('"', '\\"')

        # be sure to remove the temporary directory (if it exists) after freeview closes
        if self.tempdir and not noclean:
            command = '%s ; rm -rf %s' % (command, self.tempdir)

        # set number of OMP threads if provided
        if threads is not None:
            command = 'export OMP_NUM_THREADS=%d ; %s' % (threads, command)

        if verbose:
            print(command)

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
        extra_tags = kwargs.pop('opts', '')

        tags = ''
        for key, value in kwargs.items():

            if isinstance(value, (list, tuple)):
                value = ','.join(str(x) for x in value)

            if value is not None:
                tags += ':%s=%s' % (key, str(value))

        return tags + extra_tags

    def _vol_to_file(self, volume, name=None, force=None, ext='mgz', swap_batch_dim=False):
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

        # check if input is a numpy-convertible tensor
        if hasattr(volume, 'numpy'):
            volume = volume.numpy()

        # if input is a numpy array, convert to the appropriate fs array type
        # and let the filename creation get handled below
        if isinstance(volume, np.ndarray):

            # swap batch axis if specified
            if swap_batch_dim:
                volume = np.moveaxis(volume, 0, -1)

            volume = self._convert_ndarray(volume) if force is None else force(volume.squeeze())
            if volume is None:
                error('cannot convert array of shape %s' % str(volume.shape))
                return None

        # configure filename
        if not name:
            if isinstance(volume, Overlay):
                filename = self._unique_filename('overlay.%s' % ext)
            elif isinstance(volume, Image):
                filename = self._unique_filename('image.%s' % ext)
            else:
                filename = self._unique_filename('volume.%s' % ext)
        else:
            filename = self._unique_filename(name.replace(' ', '-') + '.' + ext)

        # ensure annotations were provided lookup tables
        if ext == 'annot' and volume.lut is None:
            error('cannot save annotation without embedded lookup table')
            return None

        # check if fs array container
        if isinstance(volume, (Overlay, Image, Volume)):
            volume.write(filename)
            return filename

        # as a final test, check if the input is possibly a nibabel image
        # we don't want nibabel to be required though, so ignore import errors
        try:
            import nibabel as nib
            if isinstance(volume, nib.spatialimages.SpatialImage):
                # convert from nib to fs volume, as it's easier to control the filename
                filename = self._unique_filename(filename)
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
            fullpath = os.path.join(directory, '%s-%02d%s' % (name, n, ext))
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
    can be existing filenames, surfaces, volumes or numpy arrays. Lists are also supported.
    Use the `Freeview` class directly to configure a more detailed session.

    Args:
        opts: Additional string of flags to add to the command.
        background: Run freeview as a background process. Defaults to True.
        swap_batch_dim: Move the first axis to the last if input is a numpy image array. Default is False.
        kwargs: kwargs are forwarded to the Freeview.show() call.
    '''
    background = kwargs.pop('background', True)
    opts = kwargs.pop('opts', '')
    swap_batch_dim = kwargs.pop('swap_batch_dim', False)

    # expand any nested lists/tuples within args
    def flatten(deep):
        for el in deep:
            if isinstance(el, (list, tuple)):
                yield from flatten(el)
            else:
                yield el

    fv = Freeview()

    for arg in flatten(args):
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
            fv.vol(arg, swap_batch_dim=swap_batch_dim)

    fv.show(background=background, opts=opts, **kwargs)


def fvoverlay(surface, overlay, background=True, opts='', verbose=False, **kwargs):
    '''Freeview wrapper to quickly load an overlay onto a surface.

    Args:
        surface: An existing surface filename.
        overlay: An existing volume filename, a numpy array, or a nibabel image to apply as an overlay.
        background: Run freeview as a background process. Defaults to True.
        verbose: Print the freeview command before running. Defaults to False.
        kwargs: kwargs are forwarded to the Freeview.show() call.
    '''
    fv = Freeview()
    fv.surf(surface, overlay=overlay)
    fv.show(background=background, opts=opts, verbose=verbose, **kwargs)
