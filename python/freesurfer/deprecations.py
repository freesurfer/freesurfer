'''
Functions and classes to be deprecated soon. Some of these routines
will raise exceptions, while others will just dispay warnings.
'''

import os
import warnings
import functools
import numpy as np
import argparse
import traceback

import freesurfer as fs



show_deps = True


def replace(replacement):
    def decorator(function):
        @functools.wraps(function)
        def wrapped(*args, **kwargs):
            
            global show_deps
            showing = show_deps and os.environ.get('FS_SURFA_PORT') is not None
            if showing:
                stack = traceback.extract_stack()[0]
                preface = ''
                if stack.filename != '<stdin>':
                    preface = fs.utils.term.dim(f'{stack.filename}  L{stack.lineno}: ')
                name = function.__name__
                if name == '__init__':
                    name = args[0].__class__.__name__
                old = fs.utils.term.blue(name)
                new = fs.utils.term.blue(replacement)
                print(f'{preface}replace {old} with {new}')
                show_deps = False

            result = function(*args, **kwargs)

            if showing:
                show_deps = True

            return result
        return wrapped
    return decorator


def deprecate(message):
    def decorator(function):
        @functools.wraps(function)
        def wrapped(*args, **kwargs):
            
            global show_deps
            showing = show_deps and os.environ.get('FS_SURFA_PORT') is not None
            if showing:
                stack = traceback.extract_stack()[0]
                preface = ''
                if stack.filename != '<stdin>':
                    preface = fs.utils.term.dim(f'{stack.filename}  L{stack.lineno}: ')
                name = function.__name__
                if name == '__init__':
                    name = args[0].__class__.__name__
                old = fs.utils.term.blue(name)
                print(f'{preface}{old} is deprecated: {message}')
                show_deps = False

            result = function(*args, **kwargs)

            if showing:
                show_deps = True

            return result
        return wrapped
    return decorator


def notneeded(function):
    @functools.wraps(function)
    def wrapped(*args, **kwargs):
        
        # export SUGGEST_SURFA_PORT=1
        global show_deps
        showing = show_deps and os.environ.get('FS_SURFA_PORT') is not None
        if showing:
            stack = traceback.extract_stack()[0]
            preface = ''
            if stack.filename != '<stdin>':
                preface = fs.utils.term.dim(f'{stack.filename}  L{stack.lineno}: ')
            name = function.__name__
            if name == '__init__':
                name = args[0].__class__.__name__
            old = fs.utils.term.blue(name)
            print(f'{preface}{old} is no longer needed')
            show_deps = False

        result = function(*args, **kwargs)

        if showing:
            show_deps = True

        return result
    return wrapped


def notimplemented(function):
    @functools.wraps(function)
    def wrapped(*args, **kwargs):
        
        # export SUGGEST_SURFA_PORT=1
        global show_deps
        showing = show_deps and os.environ.get('FS_SURFA_PORT') is not None
        if showing:
            stack = traceback.extract_stack()[0]
            preface = ''
            if stack.filename != '<stdin>':
                preface = fs.utils.term.dim(f'{stack.filename}  L{stack.lineno}: ')
            name = function.__name__
            if name == '__init__':
                name = args[0].__class__.__name__
            old = fs.utils.term.blue(name)
            print(f'{preface}{old} has no direct port to surfa: message andrew if you need this')
            show_deps = False

        result = function(*args, **kwargs)

        if showing:
            show_deps = True

        return result
    return wrapped


def unsure(function):
    @functools.wraps(function)
    def wrapped(*args, **kwargs):
        
        # export SUGGEST_SURFA_PORT=1
        global show_deps
        showing = show_deps and os.environ.get('FS_SURFA_PORT') is not None
        if showing:
            stack = traceback.extract_stack()[0]
            preface = ''
            if stack.filename != '<stdin>':
                preface = fs.utils.term.dim(f'{stack.filename}  L{stack.lineno}: ')
            name = function.__name__
            if name == '__init__':
                name = args[0].__class__.__name__
            old = fs.utils.term.blue(name)
            print(f'{preface}unsure of the best way to port {old}: message andrew if you need thigs')
            show_deps = False

        result = function(*args, **kwargs)

        if showing:
            show_deps = True

        return result
    return wrapped


# ---- parameterization ----

@replace('sf.sphere.SphericalMapBarycentric(surf).parameterize(overlay)')
def MRISP(scale, nfuncs):
    if scale == 0: scale = 1
    cols = int(round(scale * 256))
    rows = 2 * cols
    return np.zeros((cols, rows, nfuncs))

@replace('sf.sphere.SphericalMapBarycentric(surf).parameterize(overlay)')
def parameterizeSurface(surf, overlay, fno=0, mrisp=None, scale=1.0):
    return surf.parameterize(overlay)

# ---- surface geometry ----

@notneeded
def createFaceArrays(surf):
    raise DeprecationWarning('createFaceArrays is no longer - use Surface.neighboring_faces instead')

@notneeded
def initSurface(surf):
    pass

@deprecate('face normals are now auto-computed: `mesh.face_normals`')
def computeFaceNormals(surf):
    surf.compute_normals()
    return surf.face_normals

@deprecate('vertex normals are now auto-computed: `mesh.vertex_normals`')
def computeVertexNormals(surf):
    surf.compute_normals()
    return surf.vertex_normals

# ---- surface coordinate conversions ----

@unsure
def MRISsampleVertexPatch(mris, mri_normals, mri, vno, wsize):
    import nibabel as nib
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

@deprecate("just use voxsurf = surf.convert(space='vox')")
def computeVoxelCoords(surf, vol):
    return surf.surf2vox().transform(surf.vertices)

@deprecate("just use voxsurf = surf.convert(space='vox')")
def MRIStoVoxel(mris, mri):
    vox2ras = mri.get_header().get_vox2ras_tkr()
    ras2vox = npl.inv(vox2ras)
    v = np.ones((4,1))
    for vno in range(len(mris[0])):
        v[:3,:] = np.reshape(mris[0][vno], v[0:3,:].shape)
        vox = ras2vox.dot(v)
        mris[0][vno] = np.reshape(vox[:3],mris[0][vno].shape)
    return mris

@deprecate("just use surf.convert(space='vox').vertex_normals")
def MRISnormalsToVoxel(mris, mri_normals, mri):
    import nibabel as nib
    if isinstance(mri_normals, nib.freesurfer.mghformat.MGHImage):
        mri_normals = mri_normals.get_data().squeeze()
    if (len(mris) == 2):
        mris = mris[0]
    vox2ras = mri.get_header().get_vox2ras_tkr()
    ras2vox = npl.inv(vox2ras)
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

@unsure
def computeVolumeValsInNormalDirection(surf, vol, dist):
    vertex_normals = surf.vertex_normals
    voxels = surf.surf2vox().transform(surf.vertices + vertex_normals * dist)
    vi = np.rint(voxels).astype(int)
    return vol.data[vi[:,0], vi[:,1], vi[:,2], 0]

# ---- metrics ----

@deprecate('use fs.metrics.hausdorff instead')  # ATH
def hausdorffDistance(*args, **kwargs):
    return fs.metrics.hausdorff(*args, **kwargs)

# ---- utils ----

@replace('image.bbox()')
def bbox2_3D(img):
    r = np.any(img, axis=(1, 2))
    c = np.any(img, axis=(0, 2))
    z = np.any(img, axis=(0, 1))
    rmin, rmax = np.where(r)[0][[0, -1]]
    cmin, cmax = np.where(c)[0][[0, -1]]
    zmin, zmax = np.where(z)[0][[0, -1]]
    return rmin, rmax, cmin, cmax, zmin, zmax

@unsure
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

@deprecate('use argparse.ArgumentParser instead')
def ArgParser(*args, **kwargs):
    return argparse.ArgumentParser(*args, **kwargs)

def source(filename):
    raise DeprecationWarning('sourcing from python is no longer supported as it is quite dangerous')

@replace('sf.system.collect_output')
def collectOutput(command, executable='/bin/bash'):
    return fs.utils.collect_output(command, executable=executable)

@replace('sf.load_overlay(filename)')
def read_annotation(*args, **kwargs):
    pass

@notneeded
def slicing_origin(slices):
    return tuple([slc.start for slc in slices])
