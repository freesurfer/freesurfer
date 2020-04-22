"""
This module provides various import/export functions as well as the 
'levelsetsTria' function

"""

# ------------------------------------------------------------------------------

def importMGH(filename):
    """
    A function to read Freesurfer MGH files.

    Required arguments:
        - filename

    Returns:
        - vol

    Requires valid mgh file. If not found, NaNs will be returned.

    """

    import os
    import struct
    import numpy

    if not os.path.exists(filename):
        print("WARNING: could not find "+filename+", returning NaNs")
        return numpy.nan

    fp = open(filename,'rb')
    intsize = struct.calcsize('>i')
    shortsize = struct.calcsize('>h')
    floatsize = struct.calcsize('>f')
    charsize = struct.calcsize('>b')

    v = struct.unpack('>i',fp.read(intsize))[0]
    ndim1 = struct.unpack('>i',fp.read(intsize))[0]
    ndim2 = struct.unpack('>i',fp.read(intsize))[0]
    ndim3 = struct.unpack('>i',fp.read(intsize))[0]
    nframes = struct.unpack('>i',fp.read(intsize))[0]
    vtype = struct.unpack('>i',fp.read(intsize))[0]
    dof = struct.unpack('>i',fp.read(intsize))[0]

    UNUSED_SPACE_SIZE = 256
    USED_SPACE_SIZE = (3*4) + (4*3*4) # space for ras transform
    unused_space_size = UNUSED_SPACE_SIZE - 2

    ras_good_flag = struct.unpack('>h',fp.read(shortsize))[0]
    if ras_good_flag:
        # We read these in but don't process them
        # as we just want to move to the volume data
        delta = struct.unpack('>fff',fp.read(floatsize*3))
        Mdc = struct.unpack('>fffffffff',fp.read(floatsize*9))
        Pxyz_c = struct.unpack('>fff',fp.read(floatsize*3))

    unused_space_size = unused_space_size - USED_SPACE_SIZE

    for i in range(unused_space_size):
        struct.unpack('>b',fp.read(charsize))[0]

    nv = ndim1 * ndim2 * ndim3 * nframes
    vol = numpy.fromstring(fp.read(floatsize*nv),dtype=numpy.float32).byteswap()

    nvert = max([ndim1,ndim2,ndim3])
    vol = numpy.reshape(vol,(ndim1,ndim2,ndim3,nframes),order='F')
    vol = numpy.squeeze(vol)
    fp.close()

    return vol


# ------------------------------------------------------------------------------

def levelsetsTria(v, t, p, levelsets):
    """
    This is the levelsetsTria function

    """

    import numpy as np
    from scipy.sparse import csr_matrix, lil_matrix

    vLVL = list()
    lLVL = list()
    iLVL = list()

    levelsets = (np.array(levelsets, ndmin=2))

    for l in range(len(levelsets)):

        A = lil_matrix((np.shape(v)[0],np.shape(v)[0]))

        lvl = levelsets[l]

        nlvl = p[t] > lvl

        n = np.where(np.logical_or(np.sum(nlvl, axis=1) == 1 , np.sum(nlvl, axis=1) == 2))[0]

        # interpolate points

        ti = list()
        vi = list()

        for i in range(len(n)):

            # which are the outlying points in the current tria?
            oi = np.where(nlvl[n[i],:])[0]

            #  convert 2 --> 1
            if len(oi) == 2:
                oi = np.setdiff1d((0,1,2), oi)

            # find the two non - outyling points
            oix = np.setdiff1d((0, 1, 2), oi)

            # check if we have interpolated for one or both of these points before

            if np.count_nonzero(A[t[n[i], oi], t[n[i], oix[0]]]) == 0:

                # compute difference vectors between outlying point and other points

                d10 = v[ t[n[i], oix[0]], :] - v[ t[n[i], oi], :]

                # compute differences of all points to lvl to get interpolation factors

                s10 = (lvl - p[t[n[i], oi]]) / (p[t[n[i], oix[0]]] - p[t[n[i], oi]])

                # compute new points

                v10 = s10 * d10 + v[ t[n[i], oi], :]

                # update vi and index(order matters)

                vi.append(v10.tolist()[0])

                ti10 = len(vi)

                # store between which two points we are interpolating (to avoid having duplicate points)

                A[ t[n[i], oi], t[n[i], oix[0]] ] = ti10
                A[ t[n[i], oix[0]], t[n[i], oi] ] = ti10

            else:

                ti10 = int(A[ t[n[i], oi], t[n[i], oix[0]] ].toarray().item())

            # essentially the same as above, just for oix[1]

            if np.count_nonzero(A[t[n[i], oi], t[n[i], oix[1]]]) == 0:

                d20 = v[ t[n[i], oix[1]], :] - v[ t[n[i], oi], :]

                s20 = (lvl - p[t[n[i], oi]]) / (p[t[n[i], oix[1]]] - p[t[n[i], oi]])

                v20 = s20 * d20 + v[ t[n[i], oi], :]

                # update vi and index(order matters)

                vi.append(v20.tolist()[0])

                ti20 = len(vi)

                A[ t[n[i], oi], t[n[i], oix[1]] ] = ti20
                A[ t[n[i], oix[1]], t[n[i], oi] ] = ti20

            else:

                ti20 = int(A[ t[n[i], oi], t[n[i], oix[1]] ].toarray().item())

            # store new indices

            ti.append((ti10,ti20))

            # clean up

            # clear oi oix d10 d20 s10 s20 v10 v20 t10 t20

        # store

        vLVL.append(vi)
        lLVL.append(ti)
        iLVL.append(n)

    return vLVL, lLVL, iLVL
