import numpy as np
import numpy.linalg as npl
import nibabel as nib
import math

from . import sample_patch

def createFaceArrays(surf):
#    vfaces = [[]] * len(surf.vertices)
    nvertices = len(surf.vertices)
    vfaces = []
    for vno in range(nvertices):
        vfaces.append([])
    for fno, face in enumerate(surf.faces):
        for vno in face:
            vfaces[vno].append(fno)
    return vfaces

def initSurface(surf):
    if hasattr(surf,'vfaces') == False:
        surf.vfaces = createFaceArrays(surf)

def MRISP(scale, nfuncs):
    if scale == 0: scale = 1
    cols = int(round(scale * 256))
    rows = 2 * cols
    return np.zeros((cols, rows, nfuncs))


def averageRadius(surface):

    if str(type(surface)) == "<class 'freesurfer.surface.Surface'>":
        vertices = surface.vertices
    else:
        vertices = surface[0]
    xhi = yhi = zhi = -10000.0
    xlo = ylo = zlo =  10000.0

    for (x, y, z) in vertices:
        if x > xhi: xhi = x
        if x < xlo: xlo = x
        if y > yhi: yhi = y
        if y < ylo: ylo = y
        if z > zhi: zhi = z
        if z < zlo: zlo = z

    x0 = (xlo + xhi) / 2.0
    y0 = (ylo + yhi) / 2.0
    z0 = (zlo + zhi) / 2.0

    radius = 0.0
    for (x, y, z) in vertices:
        dx = x - x0;
        dy = y - y0;
        dz = z - z0;
        radius += math.sqrt(dx**2 + dy**2 + dz**2)

    radius /= len(vertices)
    return radius


def parameterizeSurface(surface, overlay, fno=0, mrisp=None, scale=1.0):

    # compute the average radius
    a = b = c = averageRadius(surface)

    if mrisp is None:
        mrisp = MRISP(scale, 1)

    mrisp.fill(0.)
    cols = mrisp.shape[0]
    rows = mrisp.shape[1]

    fillstatus = np.zeros((cols, rows))
    # enum for fillstatus
    unfilled = 0
    filling = 1
    filled = 2

    distances = np.zeros((cols, rows))
    coords = []

    # calculate total distances to a point in parameter space
    for (x, y, z) in surface[0]:
        theta = math.atan2(y/b, x/a)
        if theta < 0: theta += 2 * math.pi
        d = c**2 - z**2
        if d < 0: d = 0
        phi = math.atan2(math.sqrt(d), z)

        uf = cols * phi / math.pi
        vf = rows * theta / (2 * math.pi)
        u = int(round(uf))
        v = int(round(vf))

        if u < 0: u = -u
        if u >= cols: u = (2*cols) - (u + 1)
        if v < 0: v += rows
        if v >= rows: v -= rows

        fillstatus[u][v] = filled
        distances[u][v] += 1  # keep track of total # of nodes
        coords.append((u, v))

    # now add in curvatures proportional to their distance from the point
    for i, (u, v) in enumerate(coords):
        totalDistance = distances[u][v]
        if totalDistance > 0.0: mrisp[u][v][fno] += overlay[i] / float(totalDistance)

    # fill in values which were unmapped using soap bubble
    tmpimage = np.zeros((cols, rows, 1))
    npasses = 0
    while True:
        num_unfilled = 0
        for u in range(cols):
            for v in range(rows):
                if fillstatus[u][v] == unfilled:
                    total = 0.0
                    n = 0
                    for uk in [-1, 0, 1]:
                        u1 = u + uk
                        if u1 < 0: u1 = -u1
                        if u1 >= cols: u1 = (2*cols) - (u1 + 1)
                        for vk in [-1, 0, 1]:
                            v1 = v + vk
                            if v1 < 0: v1 += rows
                            if v1 >= rows: v1 -= rows
                            if fillstatus[u1][v1] == filled:
                                total += mrisp[u1][v1][fno]
                                n += 1
                    if n > 0:
                        total /= float(n)  # TODO float might not be necessary
                        tmpimage[u][v][fno] = total
                        fillstatus[u][v] = filling
                    else:
                        num_unfilled += 1
                else:
                    tmpimage[u][v][fno] = mrisp[u][v][fno]
        fillstatus[fillstatus == filling] = filled
        mrisp[:] = tmpimage
        npasses += 1
        if npasses > 1000:
            raise RuntimeError('could not fill parameterization')
        if num_unfilled == 0:
            break

    return mrisp


def MRIStoVoxel(mris, mri):
    vox2ras = mri.get_header().get_vox2ras_tkr()
    ras2vox = npl.inv(vox2ras)
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

def computeFaceNormals(surf):
  
  edges1 = surf.vertices[surf.faces[:,1]] - surf.vertices[surf.faces[:,0]] ;
  edges2 = surf.vertices[surf.faces[:,2]] - surf.vertices[surf.faces[:,0]] ;

  normals = np.cross(edges1,edges2)
  ln = np.linalg.norm(normals,axis=1)
  ind = np.nonzero(ln == 0)
  ln[ind] = 1
  normals /= np.transpose(np.array([ln, ln, ln]),[1,0]) 
  return(normals)
     
def computeVertexNormals(surf):
    if hasattr(surf,'vfaces') == False:
        initSurface(surf)
    face_normals = computeFaceNormals(surf)
    nvertices = surf.vertices.shape[0]
    vertex_normals = np.zeros((nvertices, 3))
    faces = surf.faces
    for vno in range(nvertices):
        ind = surf.vfaces[vno]
        vertex_normals[vno,:] = face_normals[ind[0],:].mean(axis=0)
        norm = np.linalg.norm(vertex_normals[vno,:])
        if (abs(norm) < 1e-5):
            norm = 1
        vertex_normals[vno,:] /= norm
    return vertex_normals

def computeVoxelCoords(surf, vol):
    voxels = surf.surf2vox(vol).transform(surf.vertices)
    return(voxels)

def computeVolumeValsInNormalDirection(surf, vol, dist):
    vertex_normals = computeVertexNormals(surf)
    voxels = surf.surf2vox(vol).transform(surf.vertices+vertex_normals*dist)
    vi = np.rint(voxels).astype(int)
    return vol.image[vi[:,0],vi[:,1],vi[:,2],0]

