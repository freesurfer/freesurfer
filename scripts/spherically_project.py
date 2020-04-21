#!/usr/bin/env python2
from __future__ import print_function

import optparse
import os
import sys
import nibabel.freesurfer.io as fs
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

HELPTEXT = """
Script to compute ShapeDNA using linear FEM matrices. 

After correcting sign flips, embeds a surface mesh into the spectral domain, 
then projects it onto a unit sphere.  This is scaled and rotated to match the
atlas used for FreeSurfer surface registion.


USAGE:
spherically_project  -i <input_surface> -o <output_surface>


References:

Martin Reuter et al. Discrete Laplace-Beltrami Operators for Shape Analysis and
Segmentation. Computers & Graphics 33(3):381-390, 2009
    
Martin Reuter et al. Laplace-Beltrami spectra as "Shape-DNA" of surfaces and 
solids Computer-Aided Design 38(4):342-366, 2006

Bruce Fischl at al. High-resolution inter-subject averaging and a coordinate 
system for the cortical surface. Human Brain Mapping 8:272-284, 1999

    
Dependencies:
    Python 3.5

    Scipy 0.10 or later to solve the generalized eigenvalue problem.
    http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html

    Numpy
    http://www.numpy.org

    Nibabel to read and write FreeSurfer surface meshes
    http://nipy.org/nibabel/


Original Author: Martin Reuter
Date: Jan-18-2016

Updated by Lee Tirrell to add command line options,and read/write surface files


"""

h_input = 'path to input surface'
h_output = 'path to ouput surface, spherically projected'


def options_parse():
    """
    Command line option parser for spherically_project.py
    """
    parser = optparse.OptionParser(usage=HELPTEXT)
    parser.add_option('--input',  '-i', dest='input_surf', help=h_input)
    parser.add_option('--output', '-o', dest='output_surf', help=h_output)
    (options, args) = parser.parse_args()

    if options.input_surf is None or options.output_surf is None:
        sys.exit('ERROR: Please specify input and output surfaces')

    return options


def computeAB(v, t):
    """
    computeAB(v,t) computes the two sparse symmetric matrices representing
           the Laplace Beltrami Operator for a given triangle mesh using
           the linear finite element method (assuming a closed mesh or 
           the Neumann boundary condition).

    Inputs:   v - vertices : list of lists of 3 floats
              t - triangles: list of lists of 3 int of indices (>=0) into v array

    Outputs:  A - sparse sym. (n x n) positive semi definite numpy matrix 
              B - sparse sym. (n x n) positive definite numpy matrix (inner product)

    Can be used to solve sparse generalized Eigenvalue problem: A x = lambda B x
    or to solve Poisson equation: A x = B f (where f is function on mesh vertices)
    or to solve Laplace equation: A x = 0
    or to model the operator's action on a vector x:   y = B\(Ax) 
    """  
    v = np.array(v)
    t = np.array(t)
    tnum = t.shape[0]


    # Linear local matrices on unit triangle:
    tB = (np.ones((3,3)) + np.eye(3)) / 24.0

    tA00 = np.array([[ 0.5,-0.5, 0.0],
                     [-0.5, 0.5, 0.0],
                     [ 0.0, 0.0, 0.0]])

    tA11 = np.array([[ 0.5, 0.0,-0.5],
                     [ 0.0, 0.0, 0.0],
                     [-0.5, 0.0, 0.5]])

    tA0110 = np.array([[ 1.0,-0.5,-0.5],
                       [-0.5, 0.0, 0.5],
                       [-0.5, 0.5, 0.0]])

    # Compute vertex coordinates and a difference vector for each triangle:
    v1 = v[t[:, 0], :]
    v2 = v[t[:, 1], :]
    v3 = v[t[:, 2], :]
    v2mv1 = v2 - v1
    v3mv1 = v3 - v1
    
    # Compute length^2 of v3mv1 for each triangle:
    a0 = np.sum(v3mv1 * v3mv1, axis=1)

    # Compute length^2 of v2mv1 for each triangle:
    a1 = np.sum(v2mv1 * v2mv1, axis=1)

    # Compute dot product (v2mv1*v3mv1) for each triangle:
    a0110 = np.sum(v2mv1 * v3mv1, axis=1)

    # Compute cross product and 2*vol for each triangle:
    cr  = np.cross(v2mv1,v3mv1)
    vol = np.sqrt(np.sum(cr*cr, axis=1))
    # zero vol will cause division by zero below, so set to small value:
    vol_mean = 0.001*np.mean(vol)
    vol = [vol_mean if x == 0 else x for x in vol]

    # Construct all local A and B matrices (guess: for each triangle):
    localB = np.array([x*tB for x in vol])
    localA = np.array([(1.0/xv) * (xa0*tA00 + xa1*tA11 - xa0110*tA0110) for xv, xa0, xa1, xa0110 in zip(vol,a0,a1,a0110)])

    # Construct row and col indices.
    I = np.array([np.tile(x, (3,1)) for x in t])
    J = np.array([np.transpose(np.tile(x, (3,1))) for x in t])

    # Flatten arrays and swap I and J:
    I = I.flatten()
    J = J.flatten()
    localA = localA.flatten()
    localB = localB.flatten()

    # Construct sparse matrix:
    A = sparse.csr_matrix((localA, (I, J)))
    B = sparse.csr_matrix((localB, (I, J)))

    return A, B



def laplaceTria(v, t, k=10):
    """
    Compute linear finite-element method Laplace-Beltrami spectrum
    """
    A, M = computeAB(v,t)
    eigenvalues, eigenvectors = eigsh(A, k, M, sigma=-0.01)

    #eigenvalues = eigenvalues.tolist()
    #eigenvectors = eigenvectors.tolist()
   
    return eigenvalues, eigenvectors



def sphericalProject(v, t):
    """
    spherical(v,t) computes the first three non-constant eigenfunctions
           and then projects the spectral embedding onto a sphere. This works
           when the first functions have a single closed zero level set,
           splitting the mesh into two domains each. Depending on the original 
           shape triangles could get inverted. We also flip the functions 
           according to the axes that they are aligned with for the special 
           case of brain surfaces in FreeSurfer coordiates.

    Inputs:   v - vertices : list of lists of 3 floats
              t - triangles: list of lists of 3 int of indices (>=0) into v array

    Outputs:  v - vertices : of projected mesh
              t - triangles: same as before

    """
    evecs = laplaceTria(v, t, k=4)[1]

    # flip efuncs to align to coordinates consistently
    ev1 = evecs[:,1]
    ev1maxi = np.argmax(ev1)
    ev1mini = np.argmin(ev1)
    cmax = v[ev1maxi,:]
    cmin = v[ev1mini,:]
    # axis 1 = y is aligned with this function (for brains in FS space)
    if (cmax[1] < cmin[1]):
        ev1 = -1 * ev1;
        
    ev2 = evecs[:,2]
    ev2maxi = np.argmax(ev2)
    ev2mini = np.argmin(ev2)
    cmax = v[ev2maxi,:]
    cmin = v[ev2mini,:]
    # axis 2 = z is aligned with this function (for brains in FS space)
    if (cmax[2] < cmin[2]):
        ev2 = -1 * ev2;
        
    ev3 = evecs[:,3]
    ev3maxi = np.argmax(ev3)
    ev3mini = np.argmin(ev3)
    cmax = v[ev3maxi,:]
    cmin = v[ev3mini,:]
    # axis 0 = x is aligned with this function (for brains in FS space)
    if (cmax[0] < cmin[0]):
        ev3 = -1 * ev3;

    # we map evN to -1..0..+1 (keep zero level fixed)
    # I have the feeling that this helps a little with the stretching
    # at the poles, but who knows...
    ev1min = np.amin(ev1)
    ev1max = np.amax(ev1)
    ev1[ev1<0] /= - ev1min
    ev1[ev1>0] /= ev1max

    ev2min = np.amin(ev2)
    ev2max = np.amax(ev2)
    ev2[ev2<0] /= - ev2min
    ev2[ev2>0] /= ev2max

    ev3min = np.amin(ev3)
    ev3max = np.amax(ev3)
    ev3[ev3<0] /= - ev3min
    ev3[ev3>0] /= ev3max

    # project to sphere and scaled to have the same scale/origin as FS:
    dist = np.sqrt( np.square(ev1) + np.square(ev2) + np.square(ev3))
    v[:,0] = 100 * (ev3 / dist)
    v[:,1] = 100 * (ev1 / dist)
    v[:,2] = 100 * (ev2 / dist)


    return v, t



def spherically_project_surface(insurf, outsurf):
    """ (string) -> None
    takes path to insurf, spherically projects it, outputs it to outsurf
    """
    surf = fs.read_geometry(insurf, read_metadata=True)
    projected = sphericalProject(surf[0], surf[1])
    fs.write_geometry(outsurf, projected[0], projected[1], volume_info=surf[2]) 


if __name__=="__main__":
    # Command Line options are error checking done here
    options = options_parse()
    surf_to_project = options.input_surf 
    projected_surf = options.output_surf

    print("Reading in surface: {} ...".format(surf_to_project))
    spherically_project_surface(surf_to_project, projected_surf)
    print("Outputing spherically projected surface: {}".format(projected_surf))

    sys.exit(0)


