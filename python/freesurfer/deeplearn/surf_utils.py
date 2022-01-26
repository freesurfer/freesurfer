import freesurfer as fs
import numpy as np
import sys,os
import voxelmorph_sandbox as vxms
import voxelmorph as vxm
import nibabel as nib
from  tensorflow import keras
import tensorflow as tf
from tensorflow.keras import layers as KL
from tqdm import tqdm
import copy

def build_atlas(model, wt_fname, sphere_gen, target_atlas, feature_names = ['sulc', 'curv'], pad=8):
    model.load_weights(wt_fname)
    
    features = [[] for x in range(len(feature_names))]
    sphere_gen.set_iter_no(0)
    sphere_iter = iter(sphere_gen)
    for i in range(len(sphere_gen.get_inputs())):
        input_features,input_overlay = next(sphere_iter)
        for g in range(len(feature_names)):
            outputs=model.predict([input_features,target_atlas,input_features[...,g:g+1]],verbose=0)
            feature_warped = outputs[3][0,:,:,0]
            features[g].append(feature_warped)

    features_avg = []
    features_std = []
    for g in range(len(feature_names)):
        features_avg.append(np.array(features[g]).mean(axis=0))
        features_std.append(np.array(features[g]).std(axis=0))

    atlas_means = np.transpose(np.array(features_avg), (1,2,0))[np.newaxis]
    atlas_stds = np.transpose(np.array(features_std), (1,2,0))[np.newaxis]

    return atlas_means, atlas_stds

def print_loss(model, step, training, train_loss):
    """
    Prints training progress to std. out
    :param step: iteration number
    :param training: a 0/1 indicating training/testing
    :param train_loss: model loss at current iteration
    """
    if not hasattr(print_loss, "last_step"):
        print_loss.last_step = 0  # it doesn't exist yet, so initialize it
    nsteps = step - print_loss.last_step
    print_loss.last_step = step

    s = str(step) + ": " 

    names = model.metrics_names
    if isinstance(train_loss, list) or isinstance(train_loss, np.ndarray):
        for i in range(len(train_loss)):
            s += ", " + str(train_loss[i])
    else:
        s += ", " + str(train_loss)

    s = 'step %3.3d: ' % step
    for i in range(len(train_loss)):
        if i == 3:
            continue ;  # bound loss is 0
        if i == 0:
            s += ' '
        else:
            s += ', '
        if (i > 0 and hasattr(model.outputs[i-1], '_name')):  # -1 because first entry in loss is overall loss
            name = model.outputs[i-1]._name
        else:
            name = names[i]
        s += '%s=%2.3f' % (name, train_loss[i])
        if np.isfinite(train_loss[i]) == False:
            gdb.set_trace()
    print(s)
    sys.stdout.flush()



def build_semi_model_old(vol_size, nf_enc, nf_dec, nchannels=1):
    model = network2d.miccai2018_net(vol_size,nf_enc,nf_dec, use_miccai_int=False, indexing='ij', nchannels=nchannels)
    overlay_input = keras.layers.Input(shape=vol_size + (1,))
    warp_layer = model.get_layer('diffflow').output
    warp_layer._name = 'warp'
    overlay_output = vxm.layers.SpatialTransformer(name='overlay')([overlay_input, warp_layer])
    model_semi = keras.models.Model(inputs=model.inputs+[overlay_input], outputs=model.outputs+[overlay_output])
    model_semi.outputs[0]._name = 'warped input'
    model_semi.outputs[1]._name = 'velocity'
    return model_semi

def build_semi_model(mrisp_shape, nlabels, nb_feat, nb_levels, nb_conv_per_level=2, nchannels=1):
    model = vxm.networks.VxmDense(mrisp_shape, 
                                  nb_unet_features = nb_feat, 
                                  nb_unet_levels = nb_levels, 
                                  nb_unet_conv_per_level=nb_conv_per_level, 
                                  src_feats=nchannels, 
                                  trg_feats=nchannels,
                                  int_downsize=1)
    label_input = KL.Input(shape=mrisp_shape + (nlabels,),name='label_input')
    label_output = vxm.layers.SpatialTransformer(name='warped_label')([label_input, model.references.pos_flow])
    model_semi = keras.models.Model(inputs=model.inputs+[label_input], outputs=model.outputs+[label_output])
    return model_semi

def build_semi_atlas_model(mrisp_shape, nlabels, nb_feat, nb_levels, nb_conv_per_level=2, nchannels=1):
    model = vxm.networks.VxmDense(mrisp_shape, 
                                  nb_unet_features = nb_feat, 
                                  nb_unet_levels = nb_levels, 
                                  nb_unet_conv_per_level=nb_conv_per_level, 
                                  src_feats=nchannels, 
                                  trg_feats=2*nchannels,  # means and variances of atlas
                                  int_downsize=1)
    if nlabels > 0:
        label_input = KL.Input(shape=mrisp_shape + (nlabels,),name='semi_input')
        label_output = vxm.layers.SpatialTransformer(name='warped_semi')([label_input, model.references.pos_flow])
        model_semi = keras.models.Model(inputs=model.inputs+[label_input], outputs=model.outputs+[label_output])
    else:
        model_semi = model

    return model_semi


def warp_all_inputs(model_semi, train_sphere_gen, atlas_target_mrisp_mean, pad=8):

    train_iter = iter(train_sphere_gen)
    train_sphere_gen.set_iter_no(0)
    overlays_warped = []
    features_warped = []
    for i in range(len(train_sphere_gen.get_inputs())):
        input_feature,input_overlay = next(train_iter)
        outputs=model_semi.predict([input_feature,atlas_target_mrisp_mean,input_overlay],verbose=0)
        feature_warped = outputs[0][0,pad:-pad,:,0]
        velocity = outputs[1][0,pad:-pad,:,:]
        overlay_warped = outputs[3][0,pad:-pad,:,0]

        features_warped.append(feature_warped)
        overlays_warped.append(overlay_warped)
    
    return features_warped, overlays_warped

def load_atlas(nfeats, hemi='lh', overlay_name = None, atlas_surf = None, pad=8):
    atlas_target_mrisp_mean_ov = fs.Overlay.read('%s.atlas_means.mrisp.nfeats%d.mgz' % (hemi, nfeats))
    atlas_target_mrisp_std_ov = fs.Overlay.read('%s.atlas_stds.mrisp.nfeats%d.mgz' % (hemi, nfeats))
    atlas_target_mrisp_mean = atlas_target_mrisp_mean_ov.data
    atlas_target_mrisp_std = atlas_target_mrisp_std_ov.data

    if (overlay_name is not None):
        atlas_overlay = fs.Overlay.read(overlay_name)
        atlas_overlay_mrisp = np.pad(atlas_surf.parameterize(atlas_overlay).transpose(), ((pad,pad),(0,0)),'wrap')[np.newaxis,...,np.newaxis]
    else:
        atlas_overlay_mrisp = None

    if len(atlas_target_mrisp_mean.shape) == 2:
        atlas_target_mrisp_mean = atlas_target_mrisp_mean[np.newaxis,...,np.newaxis]
    if len(atlas_target_mrisp_std.shape) == 2:
        atlas_target_mrisp_std = atlas_target_mrisp_std[np.newaxis,...,np.newaxis]
    if len(atlas_target_mrisp_mean.shape) == 3:
        atlas_target_mrisp_mean = atlas_target_mrisp_mean[np.newaxis]
    if len(atlas_target_mrisp_std.shape) == 3:
        atlas_target_mrisp_std = atlas_target_mrisp_std[np.newaxis]
    return atlas_target_mrisp_mean, atlas_target_mrisp_std, atlas_overlay_mrisp


def semi_surf_gen(mrisp_list, mrisp_label_list, batch_size=4, use_rand=True):
    nsubjects = len(mrisp_list)
    nlabels = mrisp_label_list[0].shape[-1]
    mrisp_shape = mrisp_list[0].shape[0:-1]
    nfeats = mrisp_list[0].shape[-1]
    batch_moving = np.zeros((batch_size, *mrisp_shape, nfeats))
    batch_fixed = np.zeros((batch_size, *mrisp_shape, nfeats))
    batch_moving_labels = np.zeros((batch_size, *mrisp_shape, nlabels))
    batch_fixed_labels = np.zeros((batch_size, *mrisp_shape, nlabels))
    zero_warp = np.zeros((batch_size, *mrisp_shape, 2))
    while True:
        for bno in range(batch_size):
            ind_moving = np.random.randint(0,nsubjects)
            ind_fixed = np.random.randint(0,nsubjects)
            batch_moving[bno, ...] = mrisp_list[ind_moving]
            batch_moving_labels[bno,...] = mrisp_label_list[ind_moving]

            batch_fixed[bno, ...] = mrisp_list[ind_fixed]
            batch_fixed_labels[bno,...] = mrisp_label_list[ind_fixed]

        inputs = [batch_moving, batch_fixed, batch_moving_labels]
        outputs = [batch_fixed, zero_warp, batch_fixed_labels]
        yield inputs, outputs
            
def semi_surf_template_gen(mrisp_list, mrisp_label_list, batch_size=4, use_rand=True):
    nsubjects = len(mrisp_list)
    nlabels = mrisp_label_list[0].shape[-1]
    mrisp_shape = mrisp_list[0].shape[0:-1]
    nfeats = mrisp_list[0].shape[-1]
    batch_image = np.zeros((batch_size, *mrisp_shape, nfeats))
    batch_labels = np.zeros((batch_size, *mrisp_shape, nlabels))
    zero_warp = np.zeros((batch_size, *mrisp_shape, 2))
    zero_sub_in_atlas = np.zeros((batch_size, *mrisp_shape, 1))
    ind = -1
    while True:
        for bno in range(batch_size):
            if use_rand:
                ind = np.random.randint(0,nsubjects)
            else:
                ind = np.mod(ind+1, nsubjects)

            batch_image[bno, ...] = mrisp_list[ind]
            batch_labels[bno,...] = mrisp_label_list[ind]

        inputs = [batch_image]
        outputs = [batch_image, zero_sub_in_atlas, batch_labels, zero_warp]
        yield inputs, outputs
            

def semi_surf_atlas_gen(mrisp_list, mrisp_label_list, atlas_mean, atlas_var, batch_size=4, use_rand=True):
    nsubjects = len(mrisp_list)
    nlabels = mrisp_label_list[0].shape[-1]
    mrisp_shape = mrisp_list[0].shape[0:-1]
    nfeats = mrisp_list[0].shape[-1]
    batch_moving = np.zeros((batch_size, *mrisp_shape, nfeats))
    batch_moving_labels = np.zeros((batch_size, *mrisp_shape, nlabels))
    zero_warp = np.zeros((batch_size, *mrisp_shape, 2))

    atlas_input = np.repeat(np.stack((atlas_mean, atlas_var), axis=-1)[np.newaxis], batch_size,axis=0)

    while True:
        for bno in range(batch_size):
            ind_moving = np.random.randint(0,nsubjects)
            batch_moving[bno, ...] = mrisp_list[ind_moving]
            batch_moving_labels[bno,...] = mrisp_label_list[ind_moving]


        inputs = [batch_moving, atlas_input, batch_moving_labels]
        outputs = [atlas_input, zero_warp, batch_moving_labels]
        yield inputs, outputs
            

def load_func_and_spheres(spaths, base_dir, func_dir, sdir1, sphere_name, hemi, geom_names=['sulc'], which_norm='Median', pad=0, smooth_steps=None):

    mrisps_geom = []
    mrisps_func = []
    spheres = []
    snames = []
    for spath in tqdm(spaths):
        path, sname = os.path.split(spath)
        snames.append(sname)
        sdir = os.path.join(spath, 'surf')
        fdir = os.path.join(base_dir, sdir1, sname)
        sphere_fname = os.path.join(sdir, sphere_name)
        sphere = fs.Surface.read(sphere_fname)
        geoms = []
        for geom_name in geom_names:
            geom_fname = os.path.join(sdir, hemi + '.' + geom_name)
            geom = normCurvature(nib.freesurfer.io.read_morph_data(geom_fname), which_norm='Median')
            geoms.append(geom)

        spheres.append(sphere)
        func = fs.Overlay.read(os.path.join(fdir, func_dir))
#        funcs.append(func)
        if ((sphere.get_vertex_positions().shape[0] != geoms[0].shape[0]) or
            (geoms[0].shape[0] != func.shape[0])):
            print('%s sphere (%d) does not match geom (%d) or func (%d)' % \
                  (sphere.get_vertex_positions().shape[0], geoms[0].shape[0], func.image.shape[0]))
            continue

        if smooth_steps is not None:
            func = sphere.smooth_overlay(func, smooth_steps)

        mrisp_geom = sphere.parameterize(np.transpose(np.array(geoms)))
        mrisp_geom = (mrisp_geom - mrisp_geom.mean()) / mrisp_geom.std()
        mrisp_func = sphere.parameterize(func)

        if pad > 0:
            #mrisp_geom = padSphere(mrisp_geom, pad)
            #mrisp_func = padSphere(mrisp_func, pad)
            mrisp_geom = pad_2d_image_spherically(mrisp_geom, pad_size=pad)
            mrisp_func = pad_2d_image_spherically(mrisp_func, pad_size=pad)

        if len(mrisp_func.shape) == 2:
            mrisp_func = mrisp_func[...,np.newaxis] # add a channels dimension
        if len(mrisp_geom.shape) == 2:
            mrisp_geom = mrisp_geom[...,np.newaxis] # add a channels dimension

        mrisps_geom.append(mrisp_geom)
        mrisps_func.append(mrisp_func)

    return([mrisps_geom, mrisps_func, spheres, snames])


def mrisp_semi_gen(mrisps_geom, mrisps_func, batch_size=4, use_rand=True, func_thresh=0, warp_downsize=1):
    ndata = len(mrisps_geom)
    mrisp_shape = mrisps_geom[0].shape[0:2]
    ngeom = mrisps_geom[0].shape[-1]
    nfunc = mrisps_func[0].shape[-1]
    batch_moving = np.zeros((batch_size, *mrisp_shape, ngeom))
    batch_fixed = np.zeros((batch_size, *mrisp_shape, ngeom))
    batch_moving_func = np.zeros((batch_size, *mrisp_shape, nfunc))
    batch_fixed_func = np.zeros((batch_size, *mrisp_shape, nfunc))
    zero_warp = np.zeros((batch_size, *tuple(np.array(mrisp_shape)//warp_downsize), 2))
    if use_rand == False:
        ind_moving = 0
        ind_fixed = ndata-1

    while True:
        for bno in range(batch_size):
            if use_rand:
                ind_moving = np.random.randint(0, ndata)
                ind_fixed = np.random.randint(0, ndata)
            else:
                ind_moving = np.mod(ind_moving+1, ndata)
                ind_fixed = np.mod(ind_fixed-1, ndata)

            batch_fixed[bno, ...] = mrisps_geom[ind_fixed]
            batch_moving[bno, ...] = mrisps_geom[ind_moving]
            batch_fixed_func[bno, ...] = mrisps_func[ind_fixed]
            batch_moving_func[bno, ...] = mrisps_func[ind_moving]
        
        inputs = [batch_moving, batch_fixed, batch_moving_func]
        outputs = [batch_fixed, zero_warp, batch_fixed_func]
        yield inputs, outputs

def mrisp_stacked_gen(mrisps_geom, mrisps_func, batch_size=4, use_rand=True, func_thresh=0, warp_downsize=1, bidir=False):
    ndata = len(mrisps_geom)
    mrisp_shape = mrisps_geom[0].shape[0:2]
    ngeom = mrisps_geom[0].shape[-1]
    nfunc = mrisps_geom[0].shape[-1]
    nfunc = mrisps_func[0].shape[-1]
    batch_moving = np.zeros((batch_size, *mrisp_shape, ngeom+nfunc))
    batch_fixed = np.zeros((batch_size, *mrisp_shape, ngeom+nfunc))
    zero_warp = np.zeros((batch_size, *tuple(np.array(mrisp_shape)//warp_downsize), 2))
    if use_rand == False:
        ind_moving = 0
        ind_fixed = ndata-1

    while True:
        for bno in range(batch_size):
            if use_rand:
                ind_moving = np.random.randint(0, ndata)
                ind_fixed = np.random.randint(0, ndata)
            else:
                ind_moving = np.mod(ind_moving+1, ndata)
                ind_fixed = np.mod(ind_fixed-1, ndata)

            batch_fixed[bno, ..., 0:ngeom] = mrisps_geom[ind_fixed]
            batch_fixed[bno, ..., ngeom:] = mrisps_func[ind_fixed]
            batch_moving[bno, ..., 0:ngeom] = mrisps_geom[ind_moving]
            batch_moving[bno, ..., ngeom:] = mrisps_func[ind_moving]
        
        inputs = [batch_moving, batch_fixed]
        outputs = [batch_fixed, zero_warp]
        if bidir:
            outputs += [batch_moving]
        yield inputs, outputs

def mrisp_stacked_atlas_gen(mrisps_geom, mrisps_func, mrisp_geom_mean, mrisp_func_mean, batch_size=4, use_rand=True, func_thresh=0, warp_downsize=1, bidir=False):
    ndata = len(mrisps_geom)
    mrisp_shape = mrisps_geom[0].shape[0:2]
    ngeom = mrisps_geom[0].shape[-1]
    nfunc = mrisps_func[0].shape[-1]
    batch_moving = np.zeros((batch_size, *mrisp_shape, ngeom+nfunc))
    zero_warp = np.zeros((batch_size, *tuple(np.array(mrisp_shape)//warp_downsize), 2))
    batch_atlas_geom = np.repeat(mrisp_geom_mean[np.newaxis], batch_size, axis=0)
    batch_atlas_func = np.repeat(mrisp_func_mean[np.newaxis], batch_size, axis=0)
    batch_atlas = np.zeros((batch_size, *mrisp_shape, ngeom+nfunc))
    batch_atlas[...,0:ngeom] = batch_atlas_geom
    batch_atlas[...,ngeom:] = batch_atlas_func

    if use_rand == False:
        ind_moving = 0
        ind_fixed = ndata-1

    while True:
        for bno in range(batch_size):
            if use_rand:
                ind_moving = np.random.randint(0, ndata)
            else:
                ind_moving = np.mod(ind_moving+1, ndata)

            batch_moving[bno, ..., 0:ngeom] = mrisps_geom[ind_moving]
            batch_moving[bno, ..., ngeom:] = mrisps_func[ind_moving]
        
        inputs = [batch_moving, batch_atlas]
        outputs = [batch_atlas, zero_warp]
        if bidir:
            outputs += [batch_moving]
        yield inputs, outputs

def mrisp_semi_atlas_gen(mrisps_geom, mrisps_func, mrisp_geom_mean, mrisp_func_mean, batch_size=4, use_rand=True, func_thresh=0, warp_downsize=1):
    ndata = len(mrisps_geom)
    mrisp_shape = mrisps_geom[0].shape[0:2]
    ngeom = mrisps_geom[0].shape[-1]
    nfunc = mrisps_func[0].shape[-1]

    batch_moving = np.zeros((batch_size, *mrisp_shape, ngeom))
    batch_moving_func = np.zeros((batch_size, *mrisp_shape, nfunc))

    batch_atlas_geom = np.repeat(mrisp_geom_mean[np.newaxis], batch_size, axis=0)
    batch_atlas_func = np.repeat(mrisp_func_mean[np.newaxis], batch_size, axis=0)

    zero_warp = np.zeros((batch_size, *tuple(np.array(mrisp_shape)//warp_downsize), 2))
    
    if use_rand == False:
        ind_moving = 0
        ind_fixed = ndata-1

    while True:
        for bno in range(batch_size):
            if use_rand:
                ind = np.random.randint(0, ndata)
            else:
                ind = np.mod(ind_moving+1, ndata)

            batch_moving[bno, ...] = mrisps_geom[ind]
            batch_moving_func[bno, ...] = mrisps_func[ind]
        
        inputs = [batch_moving, batch_atlas_geom, batch_moving_func]
        outputs = [batch_atlas_geom, zero_warp, batch_atlas_func]
        yield inputs, outputs


    

def parc_gen(mrisps_geom, mrisps_annot, batch_size=8, use_rand=True):
    mrisp_shape = mrisps_geom[0].shape[0:2]

    if len(mrisps_geom[0].shape) == 2:  # add feature axis
        ngeom = 1
        batch_inputs = np.zeros((batch_size, *mrisp_shape, ngeom))
    else:
        ngeom = mrisps_geom[0].shape[-1]
        batch_inputs = np.zeros((batch_size, *mrisp_shape, ngeom))

    nlabels = mrisps_annot[0].shape[-1]
    batch_outputs = np.zeros((batch_size, *mrisp_shape, nlabels))

    ind = 0
    while True:
        for bno in range(batch_size):
            if use_rand:
                ind = np.random.randint(0, len(mrisps_geom))
            else:
                ind = np.mod(ind+1, len(mrisps_geom))

            if len(mrisps_geom[0].shape) == 2:  # add feature axis
                batch_inputs[bno, ...] = mrisps_geom[ind][...,np.newaxis]
            else:
                batch_inputs[bno, ...] = mrisps_geom[ind]

            batch_outputs[bno,...] = mrisps_annot[ind]

        yield batch_inputs, batch_outputs


def fsgen(mrisps_geom, mrisp_atlas, batch_size=8, use_rand=True, warp_downsize=1, mean_stream=True):
    mrisp_shape = mrisps_geom[0].shape[0:2]

    if len(mrisps_geom[0].shape) == 2:  # add feature axis
        ngeom = 1
        batch_atlas = np.repeat(mrisp_atlas[np.newaxis,...,np.newaxis], batch_size, axis=0)
        batch_moving = np.zeros((batch_size, *mrisp_shape, ngeom))
    else:
        ngeom = mrisps_geom[0].shape[-1]
        batch_atlas = np.repeat(mrisp_atlas[np.newaxis], batch_size, axis=0)
        batch_moving = np.zeros((batch_size, *mrisp_shape, ngeom))

    zero_warp = np.zeros((batch_size, *tuple(np.array(mrisp_shape)//warp_downsize), 2))

    while True:
        for bno in range(batch_size):
            ind = np.random.randint(0, len(mrisps_geom))

            if len(mrisps_geom[0].shape) == 2:  # add feature axis
                batch_moving[bno, ...] = mrisps_geom[ind][...,np.newaxis]
            else:
                batch_moving[bno, ...] = mrisps_geom[ind]
                

        inputs = [batch_moving]
        outputs = [batch_moving, batch_atlas, zero_warp]
        if mean_stream:
            outputs += [zero_warp]

        yield inputs, outputs


def fsgen_segreg(mrisps_geom, mrisp_atlas, mrisps_annot, batch_size=8, use_rand=True, warp_downsize=1, mean_stream=True, use_logprob=False):
    mrisp_shape = mrisps_geom[0].shape[0:2]
    nclasses = mrisps_annot[0].shape[-1]

    if len(mrisps_geom[0].shape) == 2:  # add feature axis
        ngeom = 1
    else:
        ngeom = mrisps_geom[0].shape[-1]

    batch_moving = np.zeros((batch_size, *mrisp_shape, ngeom))
    batch_annot = np.zeros((batch_size, *mrisp_shape, nclasses))
    batch_atlas = np.repeat(mrisp_atlas[np.newaxis], batch_size, axis=0)

    zero_warp = np.zeros((batch_size, *tuple(np.array(mrisp_shape)//warp_downsize), 2))

    if use_logprob:  # expand annots to be [-100:100] instead of [0,1]
        mrisps_annot = copy.deepcopy(mrisps_annot)
        for mrisp in mrisps_annot:
            ind0 = np.nonzero(mrisp == 0)
            ind1 = np.nonzero(mrisp == 1)
            mrisp[ind0] = -100
            mrisp[ind1] = 0

    while True:
        for bno in range(batch_size):
            ind = np.random.randint(0, len(mrisps_geom))

            if len(mrisps_geom[0].shape) == 2:  # add feature axis
                batch_moving[bno, ...] = mrisps_geom[ind][...,np.newaxis]
            else:
                batch_moving[bno, ...] = mrisps_geom[ind]

            batch_annot[bno, ...] = mrisps_annot[ind]
                

        inputs = [batch_moving]

        if mean_stream:
            outputs = [batch_moving, batch_atlas, zero_warp, zero_warp, batch_annot]
        else:
            outputs = [batch_moving, batch_atlas, zero_warp, batch_annot]
        yield inputs, outputs



def fix_annot(annot):
    new_annot = copy.deepcopy(annot)
    nbr_annot_counts = np.zeros((int(np.ceil(annot.max()+1)),))
    for x in range(annot.shape[0]):
        for y in range(annot.shape[1]):
            nbr_annot_counts.fill(0)
            for xk in range(-1,2):
                xi = x+xk
                for yk in range(-1,2):
                    yi = y+yk
                    if xi >= annot.shape[0] or yi >= annot.shape[1]:
                        continue
                    lut_index = int(annot[xi,yi])
                    if annot[xi][yi] == lut_index:
                        nbr_annot_counts[lut_index] += 1
            old_index = int(annot[x,y])
            if not annot[x,y] == old_index or nbr_annot_counts[old_index] == 0:
                new_annot[x,y] = nbr_annot_counts.argmax()

    return new_annot.astype(np.int)
                

                                                                           

import neurite_sandbox as nes
                                                                                
def AddPatchLayers(input_model, psize=96, stride=64, nconvs=3, ksize=3):
    pos_flow =  input_model.get_layer('vxm_dense_diffflow')
    nclasses = input_model.outputs[3].get_shape().as_list()[-1]
    warped_seg = vxm.layers.SpatialTransformer(fill_value=0, interp_method='linear',name='warped_seg')([input_model.outputs[3], input_model.references.pos_flow])
    warped_image = input_model.outputs[0]
    warped_stack = KL.Concatenate(axis=-1,name='warped_stack')([warped_image, warped_seg])
    ksize = 3
    stride = 64
    psize = 96
    src_feats = warped_stack.get_shape().as_list()[-1]
    source_patches = nes.layers.Patches2D(psize=psize, stride=stride, padding='VALID', name='source_patches')(warped_stack)
    clist = []
    nfeats = warped_stack.shape[-1]
    patch_tensor = source_patches
    nconvs = 3
    for pno in range(source_patches.shape[-1]):
        prev_tensor = KL.Lambda(lambda x: x[...,pno], name='lambda_patch%d' % pno)(patch_tensor)
        for cno in range(nconvs):
            conv_output = KL.Conv2D(nfeats, (3,3), 1, padding='same', name='patch%d_conv%d' % (pno,cno), activation='relu')(prev_tensor)
            prev_tensor = conv_output

        output_shape = conv_output.get_shape().as_list()[1:]+[1]
        conv_output = KL.Reshape(output_shape)(conv_output)
        clist.append(conv_output)

    image_size = warped_stack.get_shape().as_list()[1:]
    concat = KL.Concatenate(axis=-1, name='patch_concat')(clist)
    recon_tensor = nes.layers.ImageFromPatches2D(image_size, stride=stride)(concat)
    seg_atlas = KL.Conv2D(nclasses, (3,3), 1, padding='same', name='seg_atlas', activation='linear')(recon_tensor)
    seg_image = vxm.layers.SpatialTransformer(fill_value=0, interp_method='linear', name='seg_subject')([seg_atlas, input_model.references.neg_flow])
    seg_output = KL.Activation('softmax', name='segreg_patch_out')(seg_image)
    
    print('creating final model')
    model_with_patches = tf.keras.Model(input_model.inputs, input_model.outputs[0:3]+[seg_output])
    model_linear = tf.keras.Model(model_with_patches.inputs, model_with_patches.outputs[0:3]+[seg_image])
    if hasattr(input_model, 'references'):
        model_with_patches.references = input_model.references
        model_linear.references = input_model.references
    return model_with_patches, model_linear

from tensorflow.keras import backend as K
from tensorflow.keras.layers import Layer, InputLayer, Input
import pdb as gdb

class LearnedPatches(Layer):
    """ 
    Keras Layer:  apply a set of weights that are location specific to an input
                  image. Input and output can both be multi-channel
       __init__(self, nb_channels) 
    nb_channels - # of output channels 
    only works for 2D (3D would run out of ram at the moment for reasonable sizes)
    the optional priors initialization  will init the weights based on a prior label map
    """

    def __init__(self, patch_size, nb_patches, **kwargs):
        self.nb_patches = nb_patches
        self.patch_size = patch_size

        super(LearnedPatches, self).__init__(**kwargs)

    def build(self, input_shape):
        nimage_vox = np.array(input_shape[1:]).prod()
        npatch_vox = np.array(self.patch_size).prod()
        self.kernel = self.add_weight(
            shape=(nimage_vox, npatch_vox, self.nb_patches),
            initializer='RandomNormal',
            name='LearnedPatchesKernel'
        )
#        self.bias = self.add_weight(
#            shape=(self.nb_patches,),
#            initializer='RandomNormal',
#            name='LearnedPatchesBias'
#        )
        super(LearnedPatches, self).build(input_shape)

    def get_config(self):
        config = super().get_config().copy()
        config.update({
            'nb_patches': self.nb_patches,
            'patch_size': self.patch_size,
        })
        return config

    def call(self, x):
        xlist = []
        for ch in range(self.nb_patches):
            map_fn = lambda z: tf.squeeze(tf.reshape(tf.matmul(tf.reshape(z, (1,tf.size(x))), self.kernel[...,ch]), (*self.patch_size,1)))
            xch = tf.stack(tf.map_fn(map_fn, x, dtype=tf.float32), 0)
            xlist.append(xch)

        xout = tf.stack(xlist, axis=-1)
        return xout

    def compute_output_shape(self, input_shape):
        output_shape = self.patch_size + (self.nb_patches)
        return output_shape

class LearnedUnPatches(Layer):
    """ 
    Keras Layer:  apply a set of weights that are location specific to an input
                  image. Input and output can both be multi-channel
       __init__(self, nb_channels) 
    nb_channels - # of output channels 
    only works for 2D (3D would run out of ram at the moment for reasonable sizes)
    the optional priors initialization  will init the weights based on a prior label map
    """

    def __init__(self, image_shape, **kwargs):
        self.image_shape = image_shape
        self.nimage_vox = np.array(image_shape).prod()

        super(LearnedUnPatches, self).__init__(**kwargs)

    def build(self, input_shape):
        npatch_vox = np.array(input_shape[1:-1]).prod()
        nb_patches = input_shape[-1]
        self.kernel = self.add_weight(
            shape=(npatch_vox, self.nimage_vox, nb_patches),
            initializer='RandomNormal',
            name='LearnedUnPatchesKernel'
        )
#        self.bias = self.add_weight(
#            shape=(self.nb_patches,),
#            initializer='RandomNormal',
#            name='LearnedPatchesBias'
#        )
        super(LearnedUnPatches, self).build(input_shape)

    def get_config(self):
        config = super().get_config().copy()
        config.update({
            'image_shape': self.image_shape,
        })
        return config

    def call(self, x):
        xlist = []
        for ch in range(x.shape[-1]):
            map_fn = lambda z: tf.reshape(tf.matmul(tf.reshape(z[...,ch], (1,tf.size(z[...,ch]))), self.kernel[...,ch]), self.image_shape)
            xch = tf.stack(tf.map_fn(map_fn, x, dtype=tf.float32), 0)
            xlist.append(xch)

        xout = tf.math.accumulate_n(xlist)
        return xout

    def compute_output_shape(self, input_shape):
        output_shape = self.image_shape
        return output_shape



from neurite.tf.modelio import LoadableModel, store_config_args
class PatchModel2D(LoadableModel):

    @store_config_args
    def __init__(self, inshape, nb_input_channels, nb_output_channels, nconvs=10, nfeats=15, patch_size=(16,16), nb_patches=32):
        """ 
        """
        inp = KL.Input((*inshape, nb_input_channels), name='patch_input')
        patches = LearnedPatches(patch_size, nb_patches, name='learned_patches')(inp)
        prev_tensor = patches
        for cno in range(nconvs):
            prev_tensor = KL.Conv2D(nfeats, (3,3), 1, padding='same', name='patch_conv%d' % (cno), activation='relu')(prev_tensor)

        image_recon = LearnedUnPatches((*inshape, nb_input_channels), name='unpatched')(prev_tensor)
        prev_tensor = image_recon
        for cno in range(nconvs):
            prev_tensor = KL.Conv2D(nfeats, (3,3), 1, padding='same', name='im_conv%d' % (cno), activation='relu')(prev_tensor)

        image_linear = KL.Conv2D(nb_output_channels, (3,3), 1, padding='same', name='image_linear', activation='relu')(prev_tensor)
        softmax_out = KL.Activation('softmax', name='patchout_softmax')(image_linear)

        super().__init__(inputs=[inp], outputs=[softmax_out])
        self.references = LoadableModel.ReferenceContainer()
        self.references.linear_output = image_linear


def normCurvature(curvFileName, which_norm='Median', norm_percentile=97, std_thresh=3):

    if isinstance(curvFileName, str):
        curv = fs.Overlay.read(curvFileName).data
    else:  # if not a string assume it is the curvature vector itself
        curv = curvFileName
        
    if which_norm == 'Percentile':
        normed = np.clip(curv / np.percentile(curv, norm_percentile), 0, 1)
    elif which_norm == 'Median':
        min_clip = np.percentile(curv, 100-norm_percentile)
        max_clip = np.percentile(curv, norm_percentile)
        st = np.std(np.clip(curv, min_clip, max_clip))
        normed = np.clip(((curv - np.median(curv))/st), -std_thresh, std_thresh)
    else:
        normed = (curv - curv.mean())/np.std(curv)

    return normed


def loadSphere(surfName, curvFileName, padSize=8, which_norm='Median'):

        surf = fs.Surface.read(surfName)
        curv = normCurvature(curvFileName, which_norm)

        mrisp = surf.parameterize(curv)
        cols = mrisp.shape[0]
        rows = mrisp.shape[1]

        data = mrisp.squeeze().transpose()

        #paddata = np.concatenate((data[rows-padSize:rows, :], data, data[0:padSize, :]), axis=0)

        paddata = np.pad(data, ((padSize,padSize), (0,0)), 'wrap')
        paddata = np.pad(paddata, ((0,0), (padSize,padSize)), 'reflect')
        return paddata

def padSphere(mrisp, pad):
    if len(mrisp.shape) == 2:
        paddata = np.pad(mrisp, ((pad,pad), (0,0)), 'wrap')
        paddata = np.pad(paddata, ((0,0), (pad,pad)), 'reflect')
    else:
        paddata = np.pad(mrisp, ((pad,pad), (0,0), (0,0)), 'wrap')
        paddata = np.pad(paddata, ((0,0), (pad,pad), (0,0)), 'reflect')
    return paddata

def mrisp_semi_gen(mrisps_geom, mrisps_func, batch_size=4, use_rand=True, func_thresh=0, warp_downsize=1):
    ndata = len(mrisps_geom)
    mrisp_shape = mrisps_geom[0].shape[0:2]
    ngeom = mrisps_geom[0].shape[-1]
    nfunc = mrisps_func[0].shape[-1]
    batch_moving = np.zeros((batch_size, *mrisp_shape, ngeom))
    batch_fixed = np.zeros((batch_size, *mrisp_shape, ngeom))
    batch_moving_func = np.zeros((batch_size, *mrisp_shape, nfunc))
    batch_fixed_func = np.zeros((batch_size, *mrisp_shape, nfunc))
    zero_warp = np.zeros((batch_size, *tuple(np.array(mrisp_shape)//warp_downsize), 2))
    if use_rand == False:
        ind_moving = 0
        ind_fixed = ndata-1

    while True:
        for bno in range(batch_size):
            if use_rand:
                ind_moving = np.random.randint(0, ndata)
                ind_fixed = np.random.randint(0, ndata)
            else:
                ind_moving = np.mod(ind_moving+1, ndata)
                ind_fixed = np.mod(ind_fixed-1, ndata)

            batch_fixed[bno, ...] = mrisps_geom[ind_fixed]
            batch_moving[bno, ...] = mrisps_geom[ind_moving]
            batch_fixed_func[bno, ...] = mrisps_func[ind_fixed]
            batch_moving_func[bno, ...] = mrisps_func[ind_moving]
        
        inputs = [batch_moving, batch_fixed, batch_moving_func]
        outputs = [batch_fixed, zero_warp, batch_fixed_func]
        yield inputs, outputs

def mrisp_stacked_gen(mrisps_geom, mrisps_func, batch_size=4, use_rand=True, func_thresh=0, warp_downsize=1, bidir=False):
    ndata = len(mrisps_geom)
    mrisp_shape = mrisps_geom[0].shape[0:2]
    ngeom = mrisps_geom[0].shape[-1]
    nfunc = mrisps_geom[0].shape[-1]
    nfunc = mrisps_func[0].shape[-1]
    batch_moving = np.zeros((batch_size, *mrisp_shape, ngeom+nfunc))
    batch_fixed = np.zeros((batch_size, *mrisp_shape, ngeom+nfunc))
    zero_warp = np.zeros((batch_size, *tuple(np.array(mrisp_shape)//warp_downsize), 2))
    if use_rand == False:
        ind_moving = 0
        ind_fixed = ndata-1

    while True:
        for bno in range(batch_size):
            if use_rand:
                ind_moving = np.random.randint(0, ndata)
                ind_fixed = np.random.randint(0, ndata)
            else:
                ind_moving = np.mod(ind_moving+1, ndata)
                ind_fixed = np.mod(ind_fixed-1, ndata)

            batch_fixed[bno, ..., 0:ngeom] = mrisps_geom[ind_fixed]
            batch_fixed[bno, ..., ngeom:] = mrisps_func[ind_fixed]
            batch_moving[bno, ..., 0:ngeom] = mrisps_geom[ind_moving]
            batch_moving[bno, ..., ngeom:] = mrisps_func[ind_moving]
        
        inputs = [batch_moving, batch_fixed]
        outputs = [batch_fixed, zero_warp]
        if bidir:
            outputs += [batch_moving]
        yield inputs, outputs

def mrisp_stacked_atlas_gen(mrisps_geom, mrisps_func, mrisp_geom_mean, mrisp_func_mean, batch_size=4, use_rand=True, func_thresh=0, warp_downsize=1, bidir=False):
    ndata = len(mrisps_geom)
    mrisp_shape = mrisps_geom[0].shape[0:2]
    ngeom = mrisps_geom[0].shape[-1]
    nfunc = mrisps_func[0].shape[-1]
    batch_moving = np.zeros((batch_size, *mrisp_shape, ngeom+nfunc))
    zero_warp = np.zeros((batch_size, *tuple(np.array(mrisp_shape)//warp_downsize), 2))
    batch_atlas_geom = np.repeat(mrisp_geom_mean[np.newaxis], batch_size, axis=0)
    batch_atlas_func = np.repeat(mrisp_func_mean[np.newaxis], batch_size, axis=0)
    batch_atlas = np.zeros((batch_size, *mrisp_shape, ngeom+nfunc))
    batch_atlas[...,0:ngeom] = batch_atlas_geom
    batch_atlas[...,ngeom:] = batch_atlas_func

    if use_rand == False:
        ind_moving = 0
        ind_fixed = ndata-1

    while True:
        for bno in range(batch_size):
            if use_rand:
                ind_moving = np.random.randint(0, ndata)
            else:
                ind_moving = np.mod(ind_moving+1, ndata)

            batch_moving[bno, ..., 0:ngeom] = mrisps_geom[ind_moving]
            batch_moving[bno, ..., ngeom:] = mrisps_func[ind_moving]
        
        inputs = [batch_moving, batch_atlas]
        outputs = [batch_atlas, zero_warp]
        if bidir:
            outputs += [batch_moving]
        yield inputs, outputs

def mrisp_semi_atlas_gen(mrisps_geom, mrisps_func, mrisp_geom_mean, mrisp_func_mean, batch_size=4, use_rand=True, func_thresh=0, warp_downsize=1):
    ndata = len(mrisps_geom)
    mrisp_shape = mrisps_geom[0].shape[0:2]
    ngeom = mrisps_geom[0].shape[-1]
    nfunc = mrisps_func[0].shape[-1]

    batch_moving = np.zeros((batch_size, *mrisp_shape, ngeom))
    batch_moving_func = np.zeros((batch_size, *mrisp_shape, nfunc))

    batch_atlas_geom = np.repeat(mrisp_geom_mean[np.newaxis], batch_size, axis=0)
    batch_atlas_func = np.repeat(mrisp_func_mean[np.newaxis], batch_size, axis=0)

    zero_warp = np.zeros((batch_size, *tuple(np.array(mrisp_shape)//warp_downsize), 2))
    
    if use_rand == False:
        ind_moving = 0
        ind_fixed = ndata-1

    while True:
        for bno in range(batch_size):
            if use_rand:
                ind = np.random.randint(0, ndata)
            else:
                ind = np.mod(ind_moving+1, ndata)

            batch_moving[bno, ...] = mrisps_geom[ind]
            batch_moving_func[bno, ...] = mrisps_func[ind]
        
        inputs = [batch_moving, batch_atlas_geom, batch_moving_func]
        outputs = [batch_atlas_geom, zero_warp, batch_atlas_func]
        yield inputs, outputs


    

def parc_gen(mrisps_geom, mrisps_annot, batch_size=8, use_rand=True):
    mrisp_shape = mrisps_geom[0].shape[0:2]

    if len(mrisps_geom[0].shape) == 2:  # add feature axis
        ngeom = 1
        batch_inputs = np.zeros((batch_size, *mrisp_shape, ngeom))
    else:
        ngeom = mrisps_geom[0].shape[-1]
        batch_inputs = np.zeros((batch_size, *mrisp_shape, ngeom))

    nlabels = mrisps_annot[0].shape[-1]
    batch_outputs = np.zeros((batch_size, *mrisp_shape, nlabels))

    ind = 0
    while True:
        for bno in range(batch_size):
            if use_rand:
                ind = np.random.randint(0, len(mrisps_geom))
            else:
                ind = np.mod(ind+1, len(mrisps_geom))

            if len(mrisps_geom[0].shape) == 2:  # add feature axis
                batch_inputs[bno, ...] = mrisps_geom[ind][...,np.newaxis]
            else:
                batch_inputs[bno, ...] = mrisps_geom[ind]

            batch_outputs[bno,...] = mrisps_annot[ind]

        yield batch_inputs, batch_outputs


def fsgen(mrisps_geom, mrisp_atlas, batch_size=8, use_rand=True, warp_downsize=1, mean_stream=True):
    mrisp_shape = mrisps_geom[0].shape[0:2]

    if len(mrisps_geom[0].shape) == 2:  # add feature axis
        ngeom = 1
        batch_atlas = np.repeat(mrisp_atlas[np.newaxis,...,np.newaxis], batch_size, axis=0)
        batch_moving = np.zeros((batch_size, *mrisp_shape, ngeom))
    else:
        ngeom = mrisps_geom[0].shape[-1]
        batch_atlas = np.repeat(mrisp_atlas[np.newaxis], batch_size, axis=0)
        batch_moving = np.zeros((batch_size, *mrisp_shape, ngeom))

    zero_warp = np.zeros((batch_size, *tuple(np.array(mrisp_shape)//warp_downsize), 2))

    while True:
        for bno in range(batch_size):
            ind = np.random.randint(0, len(mrisps_geom))

            if len(mrisps_geom[0].shape) == 2:  # add feature axis
                batch_moving[bno, ...] = mrisps_geom[ind][...,np.newaxis]
            else:
                batch_moving[bno, ...] = mrisps_geom[ind]
                

        inputs = [batch_moving]
        outputs = [batch_moving, batch_atlas, zero_warp]
        if mean_stream:
            outputs += [zero_warp]

        yield inputs, outputs


def fsgen_segreg(mrisps_geom, mrisp_atlas, mrisps_annot, batch_size=8, use_rand=True, warp_downsize=1, mean_stream=True, use_logprob=False):
    mrisp_shape = mrisps_geom[0].shape[0:2]
    nclasses = mrisps_annot[0].shape[-1]

    if len(mrisps_geom[0].shape) == 2:  # add feature axis
        ngeom = 1
    else:
        ngeom = mrisps_geom[0].shape[-1]

    batch_moving = np.zeros((batch_size, *mrisp_shape, ngeom))
    batch_annot = np.zeros((batch_size, *mrisp_shape, nclasses))
    batch_atlas = np.repeat(mrisp_atlas[np.newaxis], batch_size, axis=0)

    zero_warp = np.zeros((batch_size, *tuple(np.array(mrisp_shape)//warp_downsize), 2))

    if use_logprob:  # expand annots to be [-100:100] instead of [0,1]
        mrisps_annot = copy.deepcopy(mrisps_annot)
        for mrisp in mrisps_annot:
            ind0 = np.nonzero(mrisp == 0)
            ind1 = np.nonzero(mrisp == 1)
            mrisp[ind0] = -100
            mrisp[ind1] = 0

    while True:
        for bno in range(batch_size):
            ind = np.random.randint(0, len(mrisps_geom))

            if len(mrisps_geom[0].shape) == 2:  # add feature axis
                batch_moving[bno, ...] = mrisps_geom[ind][...,np.newaxis]
            else:
                batch_moving[bno, ...] = mrisps_geom[ind]

            batch_annot[bno, ...] = mrisps_annot[ind]
                

        inputs = [batch_moving]

        if mean_stream:
            outputs = [batch_moving, batch_atlas, zero_warp, zero_warp, batch_annot]
        else:
            outputs = [batch_moving, batch_atlas, zero_warp, batch_annot]
        yield inputs, outputs



def fix_annot(annot):
    new_annot = copy.deepcopy(annot)
    nbr_annot_counts = np.zeros((int(np.ceil(annot.max()+1)),))
    for x in range(annot.shape[0]):
        for y in range(annot.shape[1]):
            nbr_annot_counts.fill(0)
            for xk in range(-1,2):
                xi = x+xk
                for yk in range(-1,2):
                    yi = y+yk
                    if xi >= annot.shape[0] or yi >= annot.shape[1]:
                        continue
                    lut_index = int(annot[xi,yi])
                    if annot[xi][yi] == lut_index:
                        nbr_annot_counts[lut_index] += 1
            old_index = int(annot[x,y])
            if not annot[x,y] == old_index or nbr_annot_counts[old_index] == 0:
                new_annot[x,y] = nbr_annot_counts.argmax()

    return new_annot.astype(np.int)
                

                                                                           

import neurite_sandbox as nes
                                                                                
def AddPatchLayers(input_model, psize=96, stride=64, nconvs=3, ksize=3):
    pos_flow =  input_model.get_layer('vxm_dense_diffflow')
    nclasses = input_model.outputs[3].get_shape().as_list()[-1]
    warped_seg = vxm.layers.SpatialTransformer(fill_value=0, interp_method='linear',name='warped_seg')([input_model.outputs[3], input_model.references.pos_flow])
    warped_image = input_model.outputs[0]
    warped_stack = KL.Concatenate(axis=-1,name='warped_stack')([warped_image, warped_seg])
    ksize = 3
    stride = 64
    psize = 96
    src_feats = warped_stack.get_shape().as_list()[-1]
    source_patches = nes.layers.Patches2D(psize=psize, stride=stride, padding='VALID', name='source_patches')(warped_stack)
    clist = []
    nfeats = warped_stack.shape[-1]
    patch_tensor = source_patches
    nconvs = 3
    for pno in range(source_patches.shape[-1]):
        prev_tensor = KL.Lambda(lambda x: x[...,pno], name='lambda_patch%d' % pno)(patch_tensor)
        for cno in range(nconvs):
            conv_output = KL.Conv2D(nfeats, (3,3), 1, padding='same', name='patch%d_conv%d' % (pno,cno), activation='relu')(prev_tensor)
            prev_tensor = conv_output

        output_shape = conv_output.get_shape().as_list()[1:]+[1]
        conv_output = KL.Reshape(output_shape)(conv_output)
        clist.append(conv_output)

    image_size = warped_stack.get_shape().as_list()[1:]
    concat = KL.Concatenate(axis=-1, name='patch_concat')(clist)
    recon_tensor = nes.layers.ImageFromPatches2D(image_size, stride=stride)(concat)
    seg_atlas = KL.Conv2D(nclasses, (3,3), 1, padding='same', name='seg_atlas', activation='linear')(recon_tensor)
    seg_image = vxm.layers.SpatialTransformer(fill_value=0, interp_method='linear', name='seg_subject')([seg_atlas, input_model.references.neg_flow])
    seg_output = KL.Activation('softmax', name='segreg_patch_out')(seg_image)
    
    print('creating final model')
    model_with_patches = tf.keras.Model(input_model.inputs, input_model.outputs[0:3]+[seg_output])
    model_linear = tf.keras.Model(model_with_patches.inputs, model_with_patches.outputs[0:3]+[seg_image])
    if hasattr(input_model, 'references'):
        model_with_patches.references = input_model.references
        model_linear.references = input_model.references
    return model_with_patches, model_linear

from tensorflow.keras import backend as K
from tensorflow.keras.layers import Layer, InputLayer, Input
import pdb as gdb

class LearnedPatches(Layer):
    """ 
    Keras Layer:  apply a set of weights that are location specific to an input
                  image. Input and output can both be multi-channel
       __init__(self, nb_channels) 
    nb_channels - # of output channels 
    only works for 2D (3D would run out of ram at the moment for reasonable sizes)
    the optional priors initialization  will init the weights based on a prior label map
    """

    def __init__(self, patch_size, nb_patches, **kwargs):
        self.nb_patches = nb_patches
        self.patch_size = patch_size

        super(LearnedPatches, self).__init__(**kwargs)

    def build(self, input_shape):
        nimage_vox = np.array(input_shape[1:]).prod()
        npatch_vox = np.array(self.patch_size).prod()
        self.kernel = self.add_weight(
            shape=(nimage_vox, npatch_vox, self.nb_patches),
            initializer='RandomNormal',
            name='LearnedPatchesKernel'
        )
#        self.bias = self.add_weight(
#            shape=(self.nb_patches,),
#            initializer='RandomNormal',
#            name='LearnedPatchesBias'
#        )
        super(LearnedPatches, self).build(input_shape)

    def get_config(self):
        config = super().get_config().copy()
        config.update({
            'nb_patches': self.nb_patches,
            'patch_size': self.patch_size,
        })
        return config

    def call(self, x):
        xlist = []
        for ch in range(self.nb_patches):
            map_fn = lambda z: tf.squeeze(tf.reshape(tf.matmul(tf.reshape(z, (1,tf.size(x))), self.kernel[...,ch]), (*self.patch_size,1)))
            xch = tf.stack(tf.map_fn(map_fn, x, dtype=tf.float32), 0)
            xlist.append(xch)

        xout = tf.stack(xlist, axis=-1)
        return xout

    def compute_output_shape(self, input_shape):
        output_shape = self.patch_size + (self.nb_patches)
        return output_shape

class LearnedUnPatches(Layer):
    """ 
    Keras Layer:  apply a set of weights that are location specific to an input
                  image. Input and output can both be multi-channel
       __init__(self, nb_channels) 
    nb_channels - # of output channels 
    only works for 2D (3D would run out of ram at the moment for reasonable sizes)
    the optional priors initialization  will init the weights based on a prior label map
    """

    def __init__(self, image_shape, **kwargs):
        self.image_shape = image_shape
        self.nimage_vox = np.array(image_shape).prod()

        super(LearnedUnPatches, self).__init__(**kwargs)

    def build(self, input_shape):
        npatch_vox = np.array(input_shape[1:-1]).prod()
        nb_patches = input_shape[-1]
        self.kernel = self.add_weight(
            shape=(npatch_vox, self.nimage_vox, nb_patches),
            initializer='RandomNormal',
            name='LearnedUnPatchesKernel'
        )
#        self.bias = self.add_weight(
#            shape=(self.nb_patches,),
#            initializer='RandomNormal',
#            name='LearnedPatchesBias'
#        )
        super(LearnedUnPatches, self).build(input_shape)

    def get_config(self):
        config = super().get_config().copy()
        config.update({
            'image_shape': self.image_shape,
        })
        return config

    def call(self, x):
        xlist = []
        for ch in range(x.shape[-1]):
            map_fn = lambda z: tf.reshape(tf.matmul(tf.reshape(z[...,ch], (1,tf.size(z[...,ch]))), self.kernel[...,ch]), self.image_shape)
            xch = tf.stack(tf.map_fn(map_fn, x, dtype=tf.float32), 0)
            xlist.append(xch)

        xout = tf.math.accumulate_n(xlist)
        return xout

    def compute_output_shape(self, input_shape):
        output_shape = self.image_shape
        return output_shape



from neurite.tf.modelio import LoadableModel, store_config_args
class PatchModel2D(LoadableModel):

    @store_config_args
    def __init__(self, inshape, nb_input_channels, nb_output_channels, nconvs=10, nfeats=15, patch_size=(16,16), nb_patches=32):
        """ 
        """
        inp = KL.Input((*inshape, nb_input_channels), name='patch_input')
        patches = LearnedPatches(patch_size, nb_patches, name='learned_patches')(inp)
        prev_tensor = patches
        for cno in range(nconvs):
            prev_tensor = KL.Conv2D(nfeats, (3,3), 1, padding='same', name='patch_conv%d' % (cno), activation='relu')(prev_tensor)

        image_recon = LearnedUnPatches((*inshape, nb_input_channels), name='unpatched')(prev_tensor)
        prev_tensor = image_recon
        for cno in range(nconvs):
            prev_tensor = KL.Conv2D(nfeats, (3,3), 1, padding='same', name='im_conv%d' % (cno), activation='relu')(prev_tensor)

        image_linear = KL.Conv2D(nb_output_channels, (3,3), 1, padding='same', name='image_linear', activation='relu')(prev_tensor)
        softmax_out = KL.Activation('softmax', name='patchout_softmax')(image_linear)

        super().__init__(inputs=[inp], outputs=[softmax_out])
        self.references = LoadableModel.ReferenceContainer()
        self.references.linear_output = image_linear


def normCurvature(curvFileName, which_norm='Median', norm_percentile=97, std_thresh=3):

    if isinstance(curvFileName, str):
        curv = fs.Overlay.read(curvFileName).data
    else:  # if not a string assume it is the curvature vector itself
        curv = curvFileName
        
    if which_norm == 'Percentile':
        normed = np.clip(curv / np.percentile(curv, norm_percentile), 0, 1)
    elif which_norm == 'Median':
        min_clip = np.percentile(curv, 100-norm_percentile)
        max_clip = np.percentile(curv, norm_percentile)
        st = np.std(np.clip(curv, min_clip, max_clip))
        normed = np.clip(((curv - np.median(curv))/st), -std_thresh, std_thresh)
    else:
        normed = (curv - curv.mean())/np.std(curv)

    return normed


def loadSphere(surfName, curvFileName, padSize=8, which_norm='Median'):

        surf = fs.Surface.read(surfName)
        curv = normCurvature(curvFileName, which_norm)

        mrisp = surf.parameterize(curv)
        cols = mrisp.shape[0]
        rows = mrisp.shape[1]

        data = mrisp.squeeze().transpose()

        #paddata = np.concatenate((data[rows-padSize:rows, :], data, data[0:padSize, :]), axis=0)

        paddata = np.pad(data, ((padSize,padSize), (0,0)), 'wrap')
        paddata = np.pad(paddata, ((0,0), (padSize,padSize)), 'reflect')
        return paddata

def pad_2d_image_spherically(img, pad_size=16, mode=None, keep_batch_dim=False):
    """
    pad parameterized 2d image based on the spherical positions of its vertices
    :param img: image to pad, with a shape [batch_size, H, W, ...] or [H, W] for a single image
    :param pad_size: size of the pad
    :param mode: lat or long or None, indicating if latitude dim proceeds longitude dim or the other way around
                 if mode is None, lat or long will be inferred automatically from the shape of the img,
                 with the assumption that lat dim < long dim
    :param keep_batch_dim: binary flag indicating whether to keep the batch dimension or not for single image
    """
    if mode not in ['lat', 'long', None]:
        raise ValueError('mode should be either lat or long or None')

    # check if input is a single 2-D image or not, and pad to N-D if yes
    is_2d = False
    if img.ndim == 2:
        img = img[np.newaxis, ...]
        is_2d = True

    # set mode according to the input shape if it's None
    h, w = img.shape[1:3]
    if mode is None:
        if h < w:
            mode = 'lat'
        else:
            mode = 'long'

    # make sure lat dim is the 2nd dim and long dim is the 3rd dim
    dim_rest = list(np.arange(3, img.ndim))
    if mode == 'long':
        img = img.transpose([0, 2, 1, *dim_rest])

    if pad_size > 0:
        # pad the north pole on top
        top = img[:, 1:pad_size + 1, ...]  # get top pad without the first row (reflect)
        top = np.flip(top, axis=1)  # flip upside down
        top = np.roll(top, top.shape[2] // 2, axis=2)  # circularly shift by pi

        # similarly for the south pole on bottom
        bot = img[:, -pad_size - 1:-1, ...]
        bot = np.flip(bot, axis=1)
        bot = np.roll(bot, bot.shape[2] // 2, axis=2)

        # concatenate top and bottom before padding left and right
        img2 = np.concatenate((top, img, bot), axis=1)

        # pad left and right using the "wrap" function
        zero_pad_rest = ((0, 0),) * (img.ndim - 3)
        img3 = np.pad(img2, ((0, 0), (0, 0), (pad_size, pad_size), *zero_pad_rest), 'wrap')
    else:
        img3 = img

    # swap the dim back if it was swapped before
    if mode == 'long':
        img3 = img3.transpose([0, 2, 1, *dim_rest])

    # squeeze the tensor to 2-D if input is 2-D and the user does not want keep the extended dim
    if is_2d:
        if not keep_batch_dim:
            img3 = np.squeeze(img3)

    return img3
