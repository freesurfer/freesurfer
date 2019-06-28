import nibabel as nib
import math
import numpy as np
import freesurfer as fs
import copy

class MRISP():
    def __init__(self, surf, scale=1, fname=None):
        self.vertex_coords = None
        if (fname != None):
            mri = nib.load(fname)
            self.data = mri.get_data()
            self.affine = mri.affine
            self.scale = 256.0 / self.data.shape[0] 
        else:
            self.data = np.zeros((256.0/scale, 512.0/scale))
            self.scale = scale 
            self.affine = np.eye(4)
        self.import_surface(surf)
    def import_surface(self, surf):
        radius = fs.averageRadius(surf)
        nvertices = surf.vertices.shape[0]
        cols,rows = self.data.shape[0:2]
        self.vertex_coords = np.zeros((nvertices,2))
        vertices = copy.copy(surf.vertices)
        for vno in range(nvertices):
            x,y,z = vertices[vno,:] 
            theta = math.atan2(y/radius, x/radius)
            if theta < 0: theta += 2 * math.pi
            d = radius**2 - z**2
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
            self.vertex_coords[vno,0] = u
            self.vertex_coords[vno,1] = v
    def shape(self):
        return self.data.shape
    def get_data(self):
        assert type(self.vertex_coords) != type(None), 'surface not set in MRISP'
        nvertices = self.vertex_coords.shape[0]
        data = np.zeros((self.vertex_coords.shape[0],self.data.shape[2]))
        for vno in range(nvertices):
            u = np.rint(self.vertex_coords[vno,0]).astype(int)
            v = np.rint(self.vertex_coords[vno,1]).astype(int)
            data[vno,:] = self.data[u,v,:]
        return data



