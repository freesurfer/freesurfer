import os
import numpy as np
import math
import surfa as sf


def normCurvature(curvFileName, which_norm='Median'):
        curv = sf.load_overlay(curvFileName).data
        normed = (curv - np.median(curv))/np.std(curv)
        return normed


def spherePad(surfName, curvFileName, padSize=8, which_norm='Median'):

        surf = sf.load_mesh(surfName)
        curv = normCurvature(curvFileName, which_norm)

        mrisp = sf.sphere.SphericalMapBarycentric(surf).parameterize(curv).data
        cols = mrisp.shape[0]
        rows = mrisp.shape[1]

        data = mrisp.squeeze().transpose()

        paddata = np.concatenate((data[rows-padSize:rows, :], data, data[0:padSize, :]), axis=0)

        return paddata
