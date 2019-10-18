import os
import copy
import numpy as np
from . import bindings, LinearTransform


class Surface:

    def __init__(self, vertices, faces=None, affine=None, hemi=None):

        # TODEP - this is a temporary fix to support the previous way of loading from a file - it
        # is not an ideal way of handling things and should be removed as soon as possible
        if isinstance(vertices, str):
            print('moving foward, please load surfaces via fs.Surface.read(filename)')
            result = Surface.read(vertices)
            vertices = result.vertices
            faces = result.faces
            affine = result.affine
            hemi = result.hemi

        self.vertices = vertices
        self.faces = faces
        self.vertex_normals = None
        self.face_normals = None
        self.affine = affine
        self.hemi = hemi

        # compute face and vertex normals
        self.compute_normals()

    # ---- utility ----

    @classmethod
    def read(cls, filename):
        '''Loads surface from file.'''
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

    # ---- geometry ----

    def bbox(self):
        '''Returns the min and max coordinates that define the surface's bounding box'''
        return self.vertices.min(axis=0), self.vertices.max(axis=0)

    # ---- parameterization ----

    def parameterize(self, overlay):
        '''Parameterizes an nvertices-length array to a 256 x 512 image.
        Default parameterization method is barycentric.'''
        return bindings.surf.parameterize(self, overlay).squeeze()

    def sample_parameterization(self, image):
        '''Samples a parameterized image into an nvertices-length array.
        Default sampling method is barycentric.'''
        return bindings.surf.sample_parameterization(self, image)

    # ---- deprecations ----
    # TODEP

    def copy_geometry(self):
        '''Deprecated!'''
        raise DeprecationWarning('copy_geometry has been removed! Email andrew if you get this!!')

    def isSelfIntersecting(self):
        '''Deprecated!'''
        raise DeprecationWarning('isSelfIntersecting has been removed! Email andrew if you get this!!')

    def fillInterior(self):
        '''Deprecated!'''
        raise DeprecationWarning('fillInterior has been removed! Email andrew if you get this!!')

    def get_vertex_positions(self):
        '''Deprecated!'''
        warning('Surface.get_vertex_positions() is deprecated - use Surface.vertics directly. Sorry for the back and forth.')
        return self.vertics

    def set_vertex_positions(self, vertices):
        '''Deprecated!'''
        warning('Surface.set_vertex_positions() is deprecated - use Surface.vertics directly. Sorry for the back and forth.')
        self.vertics = vertices

    def get_face_vertices(self):
        '''Deprecated!'''
        warning('"Surface.get_face_vertices()" is deprecated - use "Surface.faces" directly. Sorry for the back and forth.')
        return self.faces

    def get_vertex_normals(self):
        '''Deprecated!'''
        warning('Surface.get_vertex_normals() is deprecated - use Surface.vertex_normals directly. Sorry for the back and forth.')
        return self.vertex_normals

    def get_vertex_faces(self):
        return ...
