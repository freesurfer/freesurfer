import numpy as np


def bbox2_3D(img):
    r = np.any(img, axis=(1, 2))
    c = np.any(img, axis=(0, 2))
    z = np.any(img, axis=(0, 1))
    rmin, rmax = np.where(r)[0][[0, -1]]
    cmin, cmax = np.where(c)[0][[0, -1]]
    zmin, zmax = np.where(z)[0][[0, -1]]
    return rmin, rmax, cmin, cmax, zmin, zmax


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
