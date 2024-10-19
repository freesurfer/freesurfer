import os
import h5py
import numpy as np
import surfa as sf
import tensorflow as tf
import voxelmorph as vxm


# Settings.
weights = {
    'joint': ('synthmorph.affine.2.h5', 'synthmorph.deform.3.h5',),
    'deform': ('synthmorph.deform.3.h5',),
    'affine': ('synthmorph.affine.2.h5',),
    'rigid': ('synthmorph.rigid.1.h5',),
}


def network_space(im, shape, center=None):
    """Construct transform from network space to the voxel space of an image.

    Constructs a coordinate transform from the space the network will operate
    in to the zero-based image index space. The network space has isotropic
    1-mm voxels, left-inferior-anterior (LIA) orientation, and no shear. It is
    centered on the field of view, or that of a reference image. This space is
    an indexed voxel space, not world space.

    Parameters
    ----------
    im : surfa.Volume
        Input image to construct the transform for.
    shape : (3,) array-like
        Spatial shape of the network space.
    center : surfa.Volume, optional
        Center the network space on the center of a reference image.

    Returns
    -------
    out : tuple of (3, 4) NumPy arrays
        Transform from network to input-image space and its inverse, thinking
        coordinates.

    """
    old = im.geom
    new = sf.ImageGeometry(
        shape=shape,
        voxsize=1,
        rotation='LIA',
        center=old.center if center is None else center.geom.center,
        shear=None,
    )

    net_to_vox = old.world2vox @ new.vox2world
    vox_to_net = new.world2vox @ old.vox2world
    return net_to_vox.matrix, vox_to_net.matrix


def transform(im, trans, shape=None, normalize=False, batch=False):
    """Apply a spatial transform to 3D image voxel data in dimensions.

    Applies a transformation matrix operating in zero-based index space or a
    displacement field to an image buffer.

    Parameters
    ----------
    im : surfa.Volume or NumPy array or TensorFlow tensor
        Input image to transform, without batch dimension.
    trans : array-like
        Transform to apply to the image. A matrix of shape (3, 4), a matrix
        of shape (4, 4), or a displacement field of shape (*space, 3),
        without batch dimension.
    shape : (3,) array-like, optional
        Output shape used for converting matrices to dense transforms. None
        means the shape of the input image will be used.
    normalize : bool, optional
        Min-max normalize the image intensities into the interval [0, 1].
    batch : bool, optional
        Prepend a singleton batch dimension to the output tensor.

    Returns
    -------
    out : float TensorFlow tensor
        Transformed image with a trailing feature dimension.

    """
    # Add singleton feature dimension if needed.
    if tf.rank(im) == 3:
        im = im[..., tf.newaxis]

    out = vxm.utils.transform(
        im, trans, fill_value=0, shift_center=False, shape=shape,
    )

    if normalize:
        out -= tf.reduce_min(out)
        out /= tf.reduce_max(out)

    if batch:
        out = out[tf.newaxis, ...]

    return out


def load_weights(model, weights):
    """Load weights into model or submodel.

    Attempts to load (all) weights into a model or one of its submodels. If
    that fails, `model` may be a submodel of what we got weights for, and we
    attempt to load the weights of a submodel (layer) into `model`.

    Parameters
    ----------
    model : TensorFlow model
        Model to initialize.
    weights : str or pathlib.Path
        Path to weights file.

    Raises
    ------
    ValueError
        If unsuccessful at loading any weights.

    """
    # Extract submodels.
    models = [model]
    i = 0
    while i < len(models):
        layers = [f for f in models[i].layers if isinstance(f, tf.keras.Model)]
        models.extend(layers)
        i += 1

    # Add models wrapping a single model in case this was done in training.
    # Requires list expansion or Python will get stuck.
    models.extend([tf.keras.Model(m.inputs, m(m.inputs)) for m in models])

    # Attempt to load all weights into one of the models.
    for mod in models:
        try:
            mod.load_weights(weights)
            return
        except ValueError as e:
            pass

    # Assume `model` is a submodel of what we got weights for.
    with h5py.File(weights, mode='r') as h5:
        layers = h5.attrs['layer_names']
        weights = [list(h5[lay].attrs['weight_names']) for lay in layers]

        # Layers with weights. Attempt loading.
        layers, weights = zip(*filter(lambda f: f[1], zip(layers, weights)))
        for lay, wei in zip(layers, weights):
            try:
                model.set_weights([h5[lay][w] for w in wei])
                return
            except ValueError as e:
                if lay is layers[-1]:
                    raise e


def register(arg):

    # Parse arguments.
    in_shape = (arg.extent,) * 3
    is_mat = arg.model in ('affine', 'rigid')

    # Threading.
    if arg.threads:
        tf.config.threading.set_inter_op_parallelism_threads(arg.threads)
        tf.config.threading.set_intra_op_parallelism_threads(arg.threads)

    # Input data.
    mov = sf.load_volume(arg.moving)
    fix = sf.load_volume(arg.fixed)
    if not len(mov.shape) == len(fix.shape) == 3:
        sf.system.fatal('input images are not single-frame volumes')

    # Transforms between native voxel and network coordinates. Voxel and network
    # spaces differ for each image. The networks expect isotropic 1-mm LIA spaces.
    # We center these on the original images, except for deformable registration:
    # this assumes prior affine registration, so we center the moving network space
    # on the fixed image, to take into account affine transforms applied by
    # resampling, updating the header, or passed on the command line alike.
    center = fix if arg.model == 'deform' else None
    net_to_mov, mov_to_net = network_space(mov, shape=in_shape, center=center)
    net_to_fix, fix_to_net = network_space(fix, shape=in_shape)

    # Coordinate transforms from and to world space. There is only one world.
    mov_to_ras = mov.geom.vox2world.matrix
    fix_to_ras = fix.geom.vox2world.matrix
    ras_to_mov = mov.geom.world2vox.matrix
    ras_to_fix = fix.geom.world2vox.matrix

    # Incorporate an initial matrix transform, mapping moving to fixed coordinates,
    # as FreeSurfer LTAs store the inverse of what you might expect. For mid-space
    # initialization, compute the matrix square root of the transform between fixed
    # and moving network space.
    if arg.init:
        init = sf.load_affine(arg.init).convert(space='voxel')
        if init.ndim != 3 \
            or not sf.transform.image_geometry_equal(mov.geom, init.source, tol=1e-3) \
            or not sf.transform.image_geometry_equal(fix.geom, init.target, tol=1e-3):
            sf.system.fatal('initial transform geometry does not match images')

        init = fix_to_net @ init @ net_to_mov
        if arg.mid_space:
            init = tf.linalg.sqrtm(init)
            if np.any(np.isnan(init)):
                sf.system.fatal(f'cannot compute matrix square root of {arg.init}')
            net_to_fix = net_to_fix @ init
            fix_to_net = np.linalg.inv(net_to_fix)

        net_to_mov = net_to_mov @ tf.linalg.inv(init)
        mov_to_net = np.linalg.inv(net_to_mov)

    # Take the input images to network space. When saving the moving image with the
    # correct voxel-to-RAS matrix after incorporating an initial matrix transform,
    # an image viewer taking this matrix into account will show an unchanged image.
    # However, the networks only see the voxel data, which have been moved.
    inputs = (
        transform(mov, net_to_mov, shape=in_shape, normalize=True, batch=True),
        transform(fix, net_to_fix, shape=in_shape, normalize=True, batch=True),
    )
    if arg.out_dir:
        os.makedirs(arg.out_dir, exist_ok=True)
        inp_1 = os.path.join(arg.out_dir, 'inp_1.mgz')
        inp_2 = os.path.join(arg.out_dir, 'inp_2.mgz')
        geom_1 = sf.ImageGeometry(in_shape, vox2world=mov_to_ras @ net_to_mov)
        geom_2 = sf.ImageGeometry(in_shape, vox2world=fix_to_ras @ net_to_fix)
        sf.Volume(inputs[0][0], geom_1).save(inp_1)
        sf.Volume(inputs[1][0], geom_2).save(inp_2)

    # Network.
    prop = dict(in_shape=in_shape, bidir=True)
    if is_mat:
        prop.update(make_dense=False, rigid=arg.model == 'rigid')
        model = vxm.networks.VxmAffineFeatureDetector(**prop)

    else:
        prop.update(mid_space=True, int_steps=arg.steps, skip_affine=arg.model == 'deform')
        model = vxm.networks.HyperVxmJoint(**prop)
        inputs = (tf.constant([arg.hyper]), *inputs)

    # Weights.
    if not arg.weights:
        fs = os.environ.get('FREESURFER_HOME')
        if not fs:
            sf.system.fatal('set environment variable FREESURFER_HOME or weights')
        arg.weights = [os.path.join(fs, 'models', f) for f in weights[arg.model]]

    for f in arg.weights:
        load_weights(model, weights=f)

    # Inference. The first transform maps from the moving to the fixed image, or
    # equivalently, from fixed to moving coordinates. The second is the inverse.
    # Convert transforms between moving and fixed network spaces to transforms
    # between the original voxel spaces.
    fw, bw = map(tf.squeeze, model(inputs))
    fw = vxm.utils.compose((net_to_mov, fw, fix_to_net), shift_center=False, shape=fix.shape)
    bw = vxm.utils.compose((net_to_fix, bw, mov_to_net), shift_center=False, shape=mov.shape)

    # Associate image geometries with the transforms. LTAs store the inverse.
    if is_mat:
        fw, bw = bw, fw
        fw = sf.Affine(fw, source=mov, target=fix, space='voxel')
        bw = sf.Affine(bw, source=fix, target=mov, space='voxel')

    else:
        fw = sf.Warp(fw, source=mov, target=fix, format=sf.Warp.Format.disp_crs)
        bw = sf.Warp(bw, source=fix, target=mov, format=sf.Warp.Format.disp_crs)

    # Output transforms.
    f = dict(space='world') if is_mat else dict(format=sf.Warp.Format.disp_ras)
    if arg.trans:
        fw.convert(**f).save(arg.trans)

    if arg.inverse:
        bw.convert(**f).save(arg.inverse)

    # Moved images.
    if arg.out_moving:
        mov.transform(fw, resample=not arg.header_only).save(arg.out_moving)

    if arg.out_fixed:
        fix.transform(bw, resample=not arg.header_only).save(arg.out_fixed)

    vmpeak = sf.system.vmpeak()
    if vmpeak is not None:
        print(f'#@# mri_synthmorph: {arg.model}, threads: {arg.threads}, VmPeak: {vmpeak}')
