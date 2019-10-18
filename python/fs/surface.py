import os
import copy
import numpy as np

from . import bindings, LinearTransform


class Surface:

    def __init__(self, vertices, faces, affine=None, hemi=None):

        # geometry
        self.vertices = vertices
        self.faces = faces
        self.affine = affine

        # face and vertex normals
        self.vertex_normals = None
        self.face_normals = None
        self.compute_normals()

        # extra parameters
        self.hemi = None

    # ---- utility ----

    @classmethod
    def read(cls, filename):
        '''Loads surface from file.'''
        return bindings.surf.read(os.path.abspath(filename))

    def write(self, filename):
        '''Writes surface to file.'''
        bindings.surf.write(self, os.path.abspath(filename))

    def copy(self):
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
        '''TODOC'''
        return bindings.surf.parameterize(self, overlay)

    def sample_parameterization(self, p):
        '''TODOC'''
        return bindings.surf.sample_parameterization(self, p)

    def parameterization_map(self):
        '''TODOC'''
        return bindings.surf.parameterization_map(self)

    # ---- deprecations ----
    # TODO remove these around mid-october (should be enough time)

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
