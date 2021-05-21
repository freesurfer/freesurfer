import os
import pdb as gdb
import nibabel as nib
import numpy as np
import freesurfer as fs
import tensorflow as tf
from tensorflow.python.eager.context import context, EAGER_MODE, GRAPH_MODE
from tensorflow.keras.callbacks import Callback
from shutil import copyfile


def switch_to_eager():
    switch_execution_mode(EAGER_MODE)


def switch_to_graph():
    switch_execution_mode(GRAPH_MODE)


def switch_execution_mode(mode):
    ctx = context()._eager_context
    ctx.mode = mode
    ctx.is_eager = mode ==  EAGER_MODE


def configure(gpu=0):
    """
    Configures the appropriate TF device from a cuda device integer.
    """
    gpuid = str(gpu)
    if gpuid is not None and (gpuid != '-1'):
        device = '/gpu:' + gpuid
        os.environ['CUDA_VISIBLE_DEVICES'] = gpuid
        # GPU memory configuration differs between TF 1 and 2
        if hasattr(tf, 'ConfigProto'):
            config = tf.ConfigProto()
            config.gpu_options.allow_growth = True
            config.allow_soft_placement = True
            tf.keras.backend.set_session(tf.Session(config=config))
        else:
            tf.config.set_soft_device_placement(True)
            for pd in tf.config.list_physical_devices('GPU'):
                tf.config.experimental.set_memory_growth(pd, True)
    else:
        device = '/cpu:0'
        os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
    
    return device



class LoopingIterator:
    def __init__(self, iterable):
        self.source = iterable
        self.iterator = iter(self.source)

    def __next__(self):
        try:
            return next(self.iterator)
        except StopIteration:
            self.iterator = iter(self.source)
            return next(self.iterator)


def pickle_generator(filenames):
    filenames = LoopingIterator(filenames)
    while True:
        yield fs.read_pickle(next(filenames))


def segment_image2D(net, im_intensity, wsize, box, nlabels=3,batch_size=256,stride=4):
    whalf = wsize/2
    im_pred = np.zeros(im_intensity.shape)
    if (box == None):
        box = [0, im_intensity.shape[0], 0, im_intensity.shape[1], 0, im_intensity.shape[2]]
    xstart = box[0]
    ystart = box[2]
    zstart = box[4]
    xend = box[1]
    yend = box[3]
    zend = box[5]
    x_test = np.zeros((batch_size, wsize, wsize, 3))
    print( 'processing z range from %d --> %d' % (zstart, zend))
    for x in np.arange(xstart, xend, stride):
        print ('x %d of %d' % (x, xend))
        for y in np.arange(ystart, yend, stride):
            i = 0
            for z in np.arange(zstart, zend):
                xmin = max(0,x-whalf)
                ymin = max(0,y-whalf)
                zmin = max(0,z-whalf)
                xmax = min(xmin+wsize, im_intensity.shape[0])
                ymax = min(ymin+wsize, im_intensity.shape[1])
                zmax = min(zmin+wsize, im_intensity.shape[2])
                xmin = xmax-wsize
                ymin = ymax-wsize
                zmin = zmax-wsize
                if (xmin<0):
                    xmax -= xmin
                    xmin = 0
                if (ymin<0):
                    ymax -= ymin
                    ymin = 0
                if (zmin<0):
                    zmax -= zmin
                    zmin = 0
                x_test[i,:,:,0] = im_intensity[x, ymin:ymax, zmin:zmax].squeeze()
                x_test[i,:,:,1] = im_intensity[xmin:xmax, y, zmin:zmax].squeeze()
                x_test[i,:,:,2] = im_intensity[xmin:xmax, ymin:ymax, z].squeeze()
                i += 1
                if (i == batch_size):
                    pred = net.predict(x_test, i, verbose=0).squeeze()
                    for xs in np.arange(stride):
                        if (x+xs >= im_pred.shape[0]):
                            break
                        for ys in np.arange(stride):
                            if (y+ys >= im_pred.shape[1]):
                                break
                            im_pred[x+xs,y+ys,z-i+1:z+1] = pred
                    i = 0
            if (i > 0):
                pred = net.predict(x_test[0:i,:], i, verbose=0).squeeze()
                for xs in np.arange(stride):
                    if (x+xs >= im_pred.shape[0]):
                        break
                    for ys in np.arange(stride):
                        if (y+ys >= im_pred.shape[1]):
                            break
                        im_pred[x+xs,y+ys,z-i+1:z+1] = pred

    return im_pred


def bbox2D(img):
    rows = np.any(img, axis=1)
    cols = np.any(img, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]

    return rmin, rmax, cmin, cmax

def bbox2_3D(img):
    r = np.any(img, axis=(1, 2))
    c = np.any(img, axis=(0, 2))
    z = np.any(img, axis=(0, 1))
    rmin, rmax = np.where(r)[0][[0, -1]]
    cmin, cmax = np.where(c)[0][[0, -1]]
    zmin, zmax = np.where(z)[0][[0, -1]]
    return rmin, rmax, cmin, cmax, zmin, zmax



def BatchGenerator3D(train_oct, train_labels, batch_size=64,wsize=32,n_labels=4, use_class_net=True, Gx=None, Gy=None,Gz=None,intensity_aug=[1.0, 1.0], augment_permute=False, unknown_whalf=None):

# save input parameters and generate circularly padded versions
    iaug_min = intensity_aug[0]
    iaug_max = intensity_aug[1]
    iaug_scale = (iaug_max-iaug_min)
    width = train_oct.shape[0]
    height = train_oct.shape[1]
    depth = train_oct.shape[2]
    whalf = int(wsize/2)
    train_oct_input = train_oct
    train_labels_input = train_labels
    im_padded = np.pad(train_oct.squeeze(), whalf, 'reflect')
    labels_padded = np.pad(train_labels.squeeze(), whalf, 'reflect')
    train_oct = im_padded
    train_labels = labels_padded

    if (use_class_net == True):
        batch_labels = np.zeros((batch_size, n_labels))
    else:
        batch_labels = np.zeros((batch_size, wsize, wsize, wsize, n_labels))
        label_patch = np.zeros((wsize,wsize,wsize,1))
    batch_intensity = np.zeros((batch_size, wsize, wsize, wsize,1))
    intensity_patch = np.zeros((wsize,wsize,wsize,1))
    found = 0

    while 1:
        if (Gx == None):
            x = np.random.randint(0, width) +whalf
            y = np.random.randint(0, height)+whalf
            z = np.random.randint(0, depth) +whalf
        else:
            x = Gx + whalf
            y = Gy + whalf
            z = Gz + whalf

        if 0:
            x = np.random.randint(180, 224) +whalf
            y = 69+whalf
            z = np.random.randint(218, 281) +whalf

        iaug = np.random.rand()*iaug_scale + iaug_min
        intensity_patch[:,:,:,0] = train_oct[x-whalf:x+whalf,y-whalf:y+whalf,z-whalf:z+whalf]*iaug
        if (unknown_whalf is not None) and (train_labels[x-unknown_whalf:x+unknown_whalf,y-unknown_whalf:y+unknown_whalf,z-unknown_whalf:z+unknown_whalf].max() < 1):
            continue  # region of all unknown near central voxel

        if augment_permute == True:
            permuted_axes = np.random.permutation((0, 1, 2))
            intensity_patch[:,:,:,0] = np.transpose(intensity_patch[:,:,:,0],permuted_axes)
        batch_intensity[found,:] = np.reshape(intensity_patch, batch_intensity[found,:].shape)
        if (use_class_net == True):
            for label in range(n_labels):
                batch_labels[found,label] = train_labels[x,y,z] == label
        else:
            label_patch[:,:,:,0] = train_labels[x-whalf:x+whalf,y-whalf:y+whalf,z-whalf:z+whalf]
            if (augment_permute == True):
                label_patch[:,:,:,0] = np.transpose(label_patch[:,:,:,0], permuted_axes)
            batch_labels[found,:] = tf.keras.utils.np_utils.to_categorical(label_patch, num_classes=n_labels) 
            
        found = found+1

        if (found >= batch_size):
            yield batch_intensity[0:found,:], batch_labels[0:found,:]
            found = 0


class WeightsSaver(Callback):
    def __init__(self, model, N, name, cp_iters = 0):
        self.model = model
        self.N = N
        self.batch = 0
        self.name = name
        self.cp_iters = cp_iters
        self.iters = 0

    def on_batch_end(self, batch, logs={}):
        if self.batch % self.N == 0:
            name = 'weights%08d.h5' % self.batch
            name = self.name
            self.model.save_weights(name)
            self.iters += 1
            if (self.cp_iters > 0 and self.iters >= self.cp_iters):
                fname,ext = os.path.splitext(name)
                cpname = fname + ".cp" + ext
                copyfile(name, cpname)
        self.batch += 1

class ModelSaver(Callback):
    def __init__(self, model, N, name, cp_iters = 0):
        self.model = model
        self.N = N
        self.batch = 0
        self.name = name
        self.cp_iters = cp_iters
        self.iters = 0

    def on_batch_end(self, batch, logs={}):
        if self.batch % self.N == 0:
#            name = '%s.%08d.h5' % (self.name,self.batch)
            name = self.name
            self.model.save(name)
            self.iters += 1
            if (self.cp_iters > 0 and self.iters >= self.cp_iters):
                fname,ext = os.path.splitext(name)
                cpname = fname + ".cp" + ext
                copyfile(name, cpname)
        self.batch += 1

def MRIStoVoxel(mris, mri):
    vox2ras = mri.get_header().get_vox2ras_tkr()
    ras2vox = np.linalg.inv(vox2ras)
    v = np.ones((4,1))
    for vno in range(len(mris[0])):
        v[:3,:] = np.reshape(mris[0][vno], v[0:3,:].shape)
        vox = ras2vox.dot(v)
        mris[0][vno] = np.reshape(vox[:3],mris[0][vno].shape)
    return mris

def MRISnormalsToVoxel(mris, mri_normals, mri):
    if isinstance(mri_normals, nib.freesurfer.mghformat.MGHImage):
        mri_normals = mri_normals.get_data().squeeze()
    if (len(mris) == 2):
        mris = mris[0]
    vox2ras = mri.get_header().get_vox2ras_tkr()
    ras2vox = np.linalg.inv(vox2ras)
    v = np.ones((4,1))
    normals = mri_normals + mris
    for vno in range(normals.shape[0]):
        v[:3,:] = np.reshape(normals[vno,:], v[0:3,:].shape)
        vox = ras2vox.dot(v)
        normals[vno,:] = np.reshape(vox[:3],normals[vno,:].shape)
        v[:3,:] = np.reshape(mris[vno,:], v[0:3,:].shape)
        vox = ras2vox.dot(v)
        normals[vno,:] -= np.reshape(vox[0:3], normals[vno,:].shape)
    return normals


def segment_surface_3D(net, mris, mri_normals, mri_intensity, wsize, batch_size=1):
    if isinstance(mri_normals, nib.freesurfer.mghformat.MGHImage):
        mri_normals = mri_normals.get_data().squeeze()
    whalf = wsize/2
    nvertices = mris.shape[0]
    found = 0
    surface_overlay = np.zeros((nvertices, 1))
    batch_intensity = np.zeros((batch_size, wsize, wsize, wsize, mri_intensity.shape[3]))
    for vno in range(nvertices):
        if (((vno+1) % 1000) == 0):
            print( 'segmenting v %d of %d: %2.2f%%' % (vno, nvertices, 100.0*vno/nvertices))
        patch = MRISsampleVertexPatch(mris, mri_normals, mri_intensity, vno, wsize)
        batch_intensity[found,:].fill(0)
        batch_intensity[found,:] = patch
        found += 1
        if (found >= batch_size):
            pred = net.predict(batch_intensity, 1, verbose=0).squeeze()
            surface_overlay[vno-found+1:vno+1] = np.reshape(pred, surface_overlay[vno-found+1:vno+1].shape)
            found = 0
    if (found > 0):
        pred = net.predict(batch_intensity[0:found,:], 1, verbose=0).squeeze()
        surface_overlay[vno-found+1:vno+1] = np.reshape(pred, surface_overlay[vno-found+1:vno+1].shape)
    return surface_overlay

def MRISsampleVertexPatch(mris, mri_normals, mri, vno, wsize):
    if isinstance(mri_normals, nib.freesurfer.mghformat.MGHImage) == True:
        mri_normals = mri_normals.get_data().squeeze()
    if (len(mris) == 2):
        mris = mris[0]
    nz = mri_normals[vno,:]
    nz /= np.linalg.norm(nz)
    if (nz[2] < 0.9):
        nx = np.cross(nz, [0,0,1])
    else:
        nx = np.cross(nz, [0,1,0])
    nx /= np.linalg.norm(nx)
    ny = np.cross(nx, nz)
    ny /= np.linalg.norm(ny)
    return sample_patch(mri, mris[vno]+(nz), nx, ny, nz, wsize)
    
def sample_patch(mri, point, nx, ny, nz, wsize):
    whalf = wsize/2
    patch = np.zeros((wsize,wsize,wsize, mri.shape[3]))
    for xk in range (-whalf,whalf):
        for yk in range (-whalf,whalf):
            for zk in range (-whalf,whalf):
                if (xk == 0 and yk == 0 and zk == 0):
                    xk = 0
                if (xk == 0 and yk == 0 and zk == whalf-1):
                    xk = 0
                xi = int(point[0] + nx[0]*xk + nx[1]*yk + nx[2]*zk)
                yi = int(point[1] + ny[0]*xk + ny[1]*yk + ny[2]*zk)
                zi = int(point[2] + nz[0]*xk + nz[1]*yk + nz[2]*zk)
                val = mri[xi,yi,zi,:]
                patch[xk+whalf,yk+whalf,zk+whalf,:] = val
    return patch
                
                
        
def segment_FCD_3D(mris, mri_normals, mri_intensity, net, wsize,batch_size=100):
    if (len(mris) == 2):
        mris = mris[0]

    nvertices = len(mris)
    pred_fcd = np.zeros((nvertices,1))
    found = 0
    for vno in range(nvertices):
        if vno and (vno % 1000)==0:
            print( 'processing vno %d of %d: %2.2f%%' % (vno, nvertices, 100.0*vno/nvertices))
        patch = MRISsampleVertexPatch(mris, mri_normals, mri_intensity, vno, wsize)
        if (vno == 0):
            patches = np.zeros((list([batch_size])+list(patch.shape)))
        patches[found,:] = patch
        found = found+1
        if (found == batch_size):
            pred = net.predict(patches, found, verbose=0)
            pred_fcd[vno-found+1:vno+1] = np.reshape(pred[:,1],(found,1))
            found = 0

    if (found > 0):
        pred = net.predict(patches, found, verbose=0)
        pred_fcd[vno-found+1:vno+1] = np.reshape(pred[0:found,1],(found,1))

    return(pred_fcd)


def segment_2D(mri_intensity, net, wsize,batch_size=1, stride=4, Gx=-1, Gy=-1, zmin=None,zmax=None):
    use_class_net = len(net.layers[-1].output_shape)<3
    n_labels = net.layers[-1].output_shape[1]
    width = mri_intensity.shape[0]
    height = mri_intensity.shape[1]
    depth = mri_intensity.shape[2]
    whalf = int(np.floor(wsize/2))
    wsize_depth = 3
    whalf_depth = wsize_depth/2
    im_padded = np.zeros((width+wsize, height+wsize, depth+2*whalf_depth))
    patches = np.zeros((batch_size, wsize, wsize, wsize_depth))
    for d in range(depth):
        im_padded[:,:,d+whalf_depth] = np.pad(mri_intensity[:,:,d].squeeze(), (whalf, whalf), 'edge')
        for d in range(whalf_depth):
            im_padded[:,:,d] = im_padded[:,:,whalf_depth]
            im_padded[:,:,depth+whalf_depth+d] = im_padded[:,:,whalf_depth+depth-1]
    pred_fcd = np.zeros((im_padded.shape[0],im_padded.shape[1], im_padded.shape[2],n_labels))

    found = 0
    if (zmax == None):
        zmax = depth+whalf_depth
    else:
        zmax += whalf_depth
    if (zmin == None):
        zmin = whalf_depth
    else:
        zmin += whalf_depth
    for z in range(zmin, zmax):
        for x in range(whalf, width+whalf,stride):
            if (x-whalf > 0 and ((x-whalf)%10)==0):
                print ('processing z %d of %d, x %d of %d' % (z-whalf_depth,zmax-whalf_depth,x-whalf, width))
            for y in range(whalf,height+whalf):
                if (x-whalf == Gx and y-whalf == Gy):
                    gdb.set_trace()
                intensity_patch = im_padded[x-whalf:x+whalf,y-whalf:y+whalf,z-whalf_depth:z+whalf_depth+1]
                
                patches[found,:] = intensity_patch
                found += 1
                if (found >= batch_size or y==height+whalf-1):
                    pred = net.predict(patches, batch_size=found,verbose=0)
                    pred_fcd[x,y-found+1:y+1,z,:] = pred[0:found,:]
                    for st in range(1,stride):
                        pred_fcd[x+st,y-found+1:y+1,z,:] = pred[0:found,:]
                    found = 0
                    
    if (use_class_net == False):
        counts[counts==0] = 1
        pred_fcd /= counts
    return_image = pred_fcd[whalf:width+whalf,whalf:height+whalf,whalf_depth:depth+whalf_depth]

    return(return_image)

def segment_unet_3D(mri_intensity, net, wsize,batch_size=1, stride=None, Gx=None, Gy=None, Gz=None, zmin=None,zmax=None):
    n_labels = net.layers[-1].output_shape[-1]
    width = mri_intensity.shape[0]
    height = mri_intensity.shape[1]
    depth = mri_intensity.shape[2]
    whalf = int(np.floor(wsize/2))
    im_padded = np.pad(mri_intensity.squeeze(), whalf, 'reflect')
    patches = np.zeros((batch_size, wsize, wsize, wsize,1))
    pred_labels = np.zeros((im_padded.shape[0],im_padded.shape[1], im_padded.shape[2],n_labels))
    wts = np.zeros((im_padded.shape[0],im_padded.shape[1], im_padded.shape[2], n_labels))
    if (stride is None):
        stride = int(wsize/4)

    intensity_patch = np.zeros((wsize,wsize,wsize,1))
    if (zmax == None):
        zmax = depth
    if (zmin == None):
        zmin = 0
    zmax += whalf
    zmin += whalf
    for x in range(whalf, width+whalf,stride):
        if (x-whalf > 0 and ((x-whalf)%(stride))==0):
            print( 'processing x %d of %d' % (x-whalf, width))
        for z in range(zmin, zmax,stride):
            found = 0
            for y in range(whalf,height+whalf, stride):
                if Gx is not None and (abs(x-whalf-Gx) < stride)  and (abs(y-whalf-Gy) < stride) and (abs(z-whalf-Gz) < stride):
                    gdb.set_trace()
                intensity_patch[:,:,:,0] = im_padded[x-whalf:x+whalf,y-whalf:y+whalf,z-whalf:z+whalf]
                patches[found,:] =  intensity_patch;
                found += 1
                if (found >= batch_size or y+stride >= height+whalf):
                    pred = net.predict(patches, batch_size=found,verbose=0)
                    for i in range(found):
                        y1 = (y - ((found-(i+1))*stride));
                        pred_labels[x-whalf:x+whalf,y1-whalf:y1+whalf,z-whalf:z+whalf,:] += pred[i,:].squeeze() ;
                        wts[x-whalf:x+whalf,y1-whalf:y1+whalf,z-whalf:z+whalf,:] += 1
                    found = 0

    wts[wts == 0] = 1
    pred_labels /= wts
    return_image = pred_labels[whalf:width+whalf,whalf:height+whalf,whalf:depth+whalf]

    return(return_image)
                

def segment_3D(mri_intensity, net, wsize,batch_size=1, stride=4, Gx=None, Gy=None, Gz=None, zmin=None,zmax=None):
    if (zmin is None):
        zmin = 0 
    if  (zmax is None):
        zmax = mri_intensity.shape[2]-1
    if (zmin >= mri_intensity.shape[2]):
        zmin = mri_intensity.shape[2]-2
    if (zmax >= mri_intensity.shape[2]):
        zmax = mri_intensity.shape[2]-1

    use_class_net = len(net.layers[-1].output_shape)<3
    if (use_class_net == False):
        return segment_unet_3D(mri_intensity, net, wsize,batch_size=batch_size, stride=stride, Gx=Gx, Gy=Gy, Gz=Gz,zmin=zmin,zmax=zmax)
    n_labels = net.layers[-1].output_shape[1]
    width = mri_intensity.shape[0]
    height = mri_intensity.shape[1]
    depth = mri_intensity.shape[2]
    whalf = int(np.floor(wsize/2))
    im_padded = np.pad(mri_intensity.squeeze(), whalf, 'reflect')
    patches = np.zeros((batch_size, wsize, wsize, wsize,1))
    pred_labels = np.zeros((im_padded.shape[0],im_padded.shape[1], im_padded.shape[2],n_labels))

    intensity_patch = np.zeros((wsize,wsize,wsize,1))
    if (zmax == None):
        zmax = depth+whalf
    else:
        zmax += whalf
    if (zmin == None):
        zmin = whalf
    else:
        zmin += whalf

    for z in range(zmin, zmax,stride):
        if ((z-whalf)%(stride))==0:
            print ('processing z %d of %d' % (z-whalf,zmax-whalf))
        for x in range(whalf, width+whalf,stride):
            if ((x-whalf)%(10*stride))==0:
                print ('processing z %d of %d, x %d of %d' % (z-whalf,zmax-whalf, x-whalf,width))
            found = 0
            for y in range(whalf,height+whalf,stride):
                if Gx is not None and (abs(x-whalf-Gx) < stride and abs(y-whalf-Gy) < stride and abs(z-whalf-Gz) < stride):
                    gdb.set_trace()
                intensity_patch[:,:,:,0] = im_padded[x-whalf:x+whalf,y-whalf:y+whalf,z-whalf:z+whalf]
                
                patches[found,:] = intensity_patch
                found += 1
                if (found >= batch_size or y+stride >= height+whalf):
                    pred = net.predict(patches, batch_size=found,verbose=0)
                    for i in range(found):
                        y1 = (y - ((found-(i+1))*stride));
                        if y1 >= pred_labels.shape[1]:
                            gdb.set_trace()
                        pred_labels[x, y1, z, :] = pred[i, :]
                    found = 0

            for y in range(whalf,height+whalf,stride):
                for stx in range(0,stride):
                    for sty in range(0,stride):
                        for stz in range(0,stride):
                            if (x-whalf == Gx and y-whalf == Gy and z-whalf == Gz):
                                gdb.set_trace()
                            pred_labels[x+stx,y+sty,z+stz,:] = pred_labels[x,y,z,:]
                    
    return_image = pred_labels[whalf:width+whalf,whalf:height+whalf,whalf:depth+whalf]

    return(return_image)

def segment_3D_lines(mri_intensity, net, wsize, nlines=1, stride=4, Gx=-1, Gy=-1, zmin=None,zmax=None):
    use_class_net = len(net.layers[-1].output_shape)<3
    n_labels = net.layers[-1].output_shape[1]
    width = mri_intensity.shape[0]
    height = mri_intensity.shape[1]
    depth = mri_intensity.shape[2]
    whalf = int(np.floor(wsize/2))
    im_padded = np.pad(mri_intensity.squeeze(), whalf, 'reflect')
    batch_size = nlines*width
    patches = np.zeros((batch_size, wsize, wsize, wsize,1))
    pred_labels = np.zeros((im_padded.shape[0],im_padded.shape[1], im_padded.shape[2],n_labels))

    intensity_patch = np.zeros((wsize,wsize,wsize,1))
    found = 0
    if (zmax == None):
        zmax = depth+whalf
    else:
        zmax += whalf
    if (zmin == None):
        zmin = whalf
    else:
        zmin += whalf
    for z in range(zmin, zmax, stride):
        if ((z-whalf)%10)==0:
            print( 'processing z %d of %d' % (z-whalf,zmax-whalf))
        lines = 0
        found = 0
        for x in range(whalf, width+whalf,stride):
            for y in range(whalf,height+whalf):
                intensity_patch[:,:,:,0] = im_padded[x-whalf:x+whalf,y-whalf:y+whalf,z-whalf:z+whalf]
                
                patches[found,:] = intensity_patch
                found += 1
        planes += 1
        if (planes >= nplanes):
            pred = net.predict(patches, batch_size=planes*width*height,verbose=0)
            pred_labels[:,:,z:z+nplanes,:] = np.reshape(pred[0:found,:], pred_labels[:,:,z:z+nplanes,:].shape)
            planes = 0
            for st in range(1,stride):
                pred_labels[:,:,z+st,:] = pred_labels[:,:,z,:]
                found = 0
        for st in range(1,stride):
            pred_labels[:,:,z+st,:] = pred_labels[:,:,z,:]
                    
    return_image = pred_labels[whalf:width+whalf,whalf:height+whalf,whalf:depth+whalf]

    return(return_image)

def bbox_3D(img):
    r = np.any(img, axis=(1, 2))
    c = np.any(img, axis=(0, 2))
    z = np.any(img, axis=(0, 1))
    rmin, rmax = np.where(r)[0][[0, -1]]
    cmin, cmax = np.where(c)[0][[0, -1]]
    zmin, zmax = np.where(z)[0][[0, -1]]
    return rmin, rmax, cmin, cmax, zmin, zmax

def bbox2(img):
    rows = np.any(img, axis=1)
    cols = np.any(img, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]

    return rmin, rmax, cmin, cmax




def sort_batch(input_batch, input_labels, probe_patch):
    sq = (input_batch-probe_patch)**2
    rms = np.sqrt(np.average(sq,axis=(1,2,3)))
    ind = np.argsort(rms)
    return input_batch[ind,:], input_labels[ind,:]

def compute_intensity_stats(paths, vname_list, method='mean', ranges=None, vol_list = None):
    nvols = len(vname_list)
    means = np.zeros(nvols)
    stds = np.zeros(nvols)
    for pno, path in enumerate(paths):
        for vno, vname in enumerate(vname_list):
            fname = os.path.join(path, vname)
            if (vol_list == None):
                print('%d of %d: file %s' % (pno, len(paths), fname))
                mri = fs.Volume(fname)
            else:
                mri = vol_list[pno][vno]
            im = mri.image.astype('float64')
            if (method == 'mean'):
                val = im.mean()
            elif method == 'histo':
                if ranges is not None:
                    ind = np.where(np.logical_and(im > ranges[vno,0], im < ranges[vno,1]))  
                    im = im[ind]
                try: 
                    hist,edges = np.histogram(im.flatten(),bins='auto')
                except ValueError:
                    print('exception')
                    hist,edges = np.histogram(im.flatten(),bins=500)
                    
                val = edges[hist.argmax()]
                del edges, hist
            shape = mri.image.shape
            if (vol_list == None):
                del mri
            means[vno] += val
            stds[vno] += (val*val)
    
    for vno in range(nvols):
        means[vno] /= len(paths)
        stds[vno] = np.sqrt(stds[vno]/len(paths) - means[vno]*means[vno])

    return means, stds, shape

def histo_norm_intensities_to_mode(vol, wsize = None,nbins=100, target_val=1.0):
# place the mode of the histo at the target val
    if wsize is None:
        wsize = int(min(vol.shape[0:3])/2)
    w,h,d = vol.shape[0:3]
    whalf = int(wsize/2)
    x0=int(max(0,w/2-whalf))
    y0=int(max(0,h/2-whalf))
    z0=int(max(0,d/2-whalf))
    x1 = int(min(w,w/2+wsize))
    y1 = int(min(w,h/2+wsize))
    z1 = int(min(w,d/2+wsize))
    histo = np.histogram(vol[x0:x1,y0:y1,z0:z1],bins=nbins)
    if (histo[1][0] == 0):
        histo[0][0] = 0   # remove 0 bin if present
    max_bin = np.argwhere(histo[0] == max(histo[0]))[0][0]
    max_val = histo[1][max_bin]
    if max_val <= 0:
        max_val = 1
    return vol * (target_val / max_val)

def histo_norm_intensities(vol,nbins=100, target_range=[0.0, 1.0], anchors=[0.1,0.9], force_zero_to_zero=True):
# rescale intensities to 90th percentile got to 90th % of range and 10th to 10th
# if force_zero_to_zero is true ignore bottom end of range and force 0 to map to 0
    histo = np.histogram(vol,bins=nbins)
    if (histo[1][0] == 0):
        histo[0][0] = 0   # remove 0 bin if present
    cumhisto = np.cumsum(histo[0]/histo[0].sum()) # cumulative as % of total (cdf)
    histo_bins = histo[1][0:len(histo[0])]
    high_bin = -1
    low_bin = -1
    for bin,val in enumerate(histo_bins):
        if low_bin < 0 and cumhisto[bin] > anchors[0]:
            low_bin = bin
        if high_bin < 0 and cumhisto[bin] > anchors[1]:
            hi_bin = bin
    trange = target_range[1]-target_range[0]
    if (force_zero_to_zero):
        x1 = 0
        y1 = 0
    else:
        x1 = histo_bins[low_bin]
        y1 = trange * anchors[0] + target_range[0]
    y2 = trange * anchors[0] + target_range[1]
    x2 = histo_bins[high_bin]
    m = (y2-y1) / (x2-x1)
    b = y1 - m*x1
            
    return vol * m + b


  
def mat_print(A):
    if A.ndim==1:
        print(A)
    else:
        w = max([len(str(s)) for s in A]) 
        print(u'\u250c'+u'\u2500'*w+u'\u2510') 
        for AA in A:
            print(' ', end='')
            print('[', end='')
            for i,AAA in enumerate(AA[:-1]):
                w1=max([len(str(s)) for s in A[:,i]])
                print(str(AAA)+' '*(w1-len(str(AAA))+1),end='')
            w1=max([len(str(s)) for s in A[:,-1]])
            print(str(AA[-1])+' '*(w1-len(str(AA[-1]))),end='')
            print(']')
        print(u'\u2514'+u'\u2500'*w+u'\u2518')  


def set_trainable(model, trainable):
    model.trainable=trainable
    for l in model.layers:
        if hasattr(l, 'layers'):
            set_trainable(l, trainable)
        l.trainable = trainable
    

def rebase_labels(labels):
    '''Rebase labels and return lookup table (LUT) to convert to new labels in
    interval [0, N[ as: LUT[label_map]. Be sure to pass all possible labels.'''
    labels = np.unique(labels) # Sorted.
    assert np.issubdtype(labels.dtype, np.integer), 'non-integer data'
    lab_to_ind = np.zeros(np.max(labels) + 1, dtype='int_')
    for i, lab in enumerate(labels):
        lab_to_ind[ lab ] = i
    ind_to_lab = labels
    return lab_to_ind, ind_to_lab

def np_one_hot(targets, nb_classes):
    res = np.eye(nb_classes)[np.array(targets.astype(np.int)).reshape(-1)]
    return res.reshape(list(targets.shape)+[nb_classes])



