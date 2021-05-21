import os
import copy
import numpy as np

from . import bindings, warning, Image, Overlay
from .transform import Transformable, LinearTransform


class Surface(Transformable):
    '''Triangular mesh topology represented by arrays of vertices and faces.'''

    def __init__(self, vertices, faces=None, hemi=None, geom=None):

        # make sure a string isn't being provided
        if isinstance(vertices, str):
            raise ValueError('if loading from file, use the Surface.read() class method')

        self.vertices = vertices
        self.faces = faces
        self.vertex_normals = None
        self.face_normals = None
        self.vertex_tangents_1 = None
        self.vertex_tangents_2 = None
        self.hemi = hemi
        self.geom = geom

        # compute face and vertex normals
        self.compute_normals()

    def __eq__(self, other):
        '''Check for equality.'''
        equal = np.array_equal(self.vertices, other.vertices) and \
                np.array_equal(self.faces, other.faces) and \
                self.hemi == other.hemi and \
                self.geom == other.geom
        return equal

    # ---- utility ----

    @classmethod
    def read(cls, filename):
        '''Loads surface from file.'''
        if not os.path.isfile(filename):
            raise ValueError('surface file %s does not exist' % filename)
        return bindings.surf.read(os.path.abspath(filename))

    def write(self, filename):
        '''Writes surface to file.'''
        bindings.surf.write(self, os.path.abspath(filename))

    def copy(self):
        '''Returns a deep copy of the surface.'''
        return copy.deepcopy(self)

    # ---- mesh topology ----

    @property
    def nvertices(self):
        '''Total number of vertices in the mesh.'''
        return len(self.vertices)

    @property
    def nfaces(self):
        '''Total number of faces in the mesh.'''
        return len(self.faces)

    def compute_normals(self):
        '''Compute and cache face and vertex normals.'''
        bindings.surf.compute_normals(self)
        # reset vertex tangents since normals have been updates
        self.vertex_tangents_1 = None
        self.vertex_tangents_2 = None

    def compute_tangents(self):
        '''Compute and cache vertex tangents along primary curvature direction.'''
        bindings.surf.compute_tangents(self)

    def compute_euler(self):
        '''Computes euler number of the mesh.'''
        return bindings.surf.compute_euler(self)

    def neighboring_faces(self, vertex):
        '''List of face indices that neighbor a vertex.'''
        return np.where(self.faces == vertex)[0]

    def neighboring_vertices(self, vertex):
        '''List of vertices that immediately neighbor a vertex.'''
        neighbors = np.unique([self.faces[fno] for fno in self.neighboring_faces(vertex)])
        return neighbors[neighbors != vertex]  # be sure to remove the central vertex

    # ---- geometry ----

    def bbox(self):
        '''Returns the (min, max) coordinates that define the surface's bounding box.'''
        return self.vertices.min(axis=0), self.vertices.max(axis=0)

    def geometry(self):
        '''Returns the geometry associated with the source volume.'''
        return self.geom

    def transform(self, lt):
        '''
        Returns a realigned surface in a coordinate space defined by
        the provided LinearTransform. For example, to convert surface
        vertices from fs surface coordinates to voxel coordinates:

        xform = surf.surf2vox()
        voxelsurf = surf.transform(xform)
        '''
        vertices = LinearTransform.ensure(lt).transform(self.vertices)
        return Surface(vertices, self.faces, hemi=self.hemi, geom=self.geom)

    def surf2vox(self, vol=None):
        '''
        LinearTransform that maps surface coordinates to crs coordinates.
        Target coordinates correspond to the source volume geometry, unless
        the optional vol argument is provided, in which case it will map
        to a different volume.
        '''
        if vol is not None:
            return LinearTransform.matmul(vol.geometry().ras2vox(), self.surf2ras())
        else:
            return super().surf2vox()

    def vox2surf(self, vol=None):
        '''
        LinearTransform that maps crs coordinates to surface coordinates.
        Source coordinates are assumed to be in source volume space, unless
        the optional vol argument is provided, in which case it will map
        from a different volume.
        '''
        if vol is not None:
            return LinearTransform.matmul(self.ras2surf(), vol.geometry().vox2ras())
        else:
            return super().vox2surf()

    # ---- overlay utils ----

    def parameterize(self, overlay, scale=1, interp='barycentric'):
        '''
        Parameterizes an nvertices-length overlay to an image. Interpolation method can be
        'barycentric' (default) or 'nearest'.
        '''
        interp = interp.lower()
        if interp not in ('barycentric', 'nearest'):
            raise ValueError('%s is not a valid interpolation method' % interp)
        
        data = Overlay.ensure(overlay).data
        if len(data) != self.nvertices:
            raise ValueError('overlay length (%d) differs from vertex count (%d)' % (len(data), self.nvertices))
        
        param = bindings.surf.parameterize(self, data, scale, interp).squeeze()
        
        if interp == 'nearest':
            param = param.astype(data.dtype)
        return param

    def sample_parameterization(self, image, interp='barycentric'):
        '''
        Samples a parameterized image into an nvertices-length array. Sampling method can be
        'barycentric' (default) or 'nearest'.
        '''
        interp = interp.lower()
        if interp not in ('barycentric', 'nearest'):
            raise ValueError('%s is not a valid interpolation method' % interp)

        data = Image.ensure(image).data
        overlay = bindings.surf.sample_parameterization(self, data, interp.lower())

        if interp == 'nearest':
            overlay = overlay.astype(data.dtype)
        return overlay

    def smooth_overlay(self, overlay, steps):
        '''Smooths an overlay along the mesh vertices.'''
        overlay = Overlay.ensure(overlay)
        return bindings.surf.smooth_overlay(self, overlay, steps)

    # ---- deprecations ----

    def parameterization_map(self):  # TODEP
        '''Deprecated! '''
        raise DeprecationWarning('parameterization_map() has been removed! Email andrew!') 

    def copy_geometry(self, vol):  # TODEP
        '''Deprecated! '''
        warning('copy_geometry(vol) has been removed! Just use: surf.geom = vol.geometry()')
        self.geom = vol.geometry()

    def isSelfIntersecting(self):  # TODEP
        '''Deprecated!'''
        raise DeprecationWarning('isSelfIntersecting has been removed! Email andrew if you get this!!')

    def fillInterior(self):  # TODEP
        '''Deprecated!'''
        raise DeprecationWarning('fillInterior has been removed! Email andrew if you get this!!')

    def get_vertex_positions(self):  # TODEP
        '''Deprecated - use Surface.vertices directly'''
        return self.vertices

    def set_vertex_positions(self, vertices):  # TODEP
        '''Deprecated - use Surface.vertices directly'''
        self.vertices = vertices

    def get_face_vertices(self):  # TODEP
        '''Deprecated - use Surface.faces directly'''
        warning('"Surface.get_face_vertices()" is deprecated - use "Surface.faces" directly. Sorry for the back and forth...')
        return self.faces

    def get_vertex_normals(self):  # TODEP
        '''Deprecated - use Surface.vertex_normals directly'''
        warning('"Surface.get_vertex_normals()"" is deprecated - use "Surface.vertex_normals" directly. Sorry for the back and forth...')
        return self.vertex_normals

    def get_vertex_faces(self):  # TODEP
        '''Deprecated - use Surface.neighboring_faces instead'''
        raise DeprecationWarning('get_vertex_faces has been removed! Use Surface.neighboring_faces or email andrew if you get this!!!!')
