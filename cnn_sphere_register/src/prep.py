import os
import numpy as np
import math

import freesurfer as fs


def normCurvature(curvFileName, which_norm='Median'):

        curv = fs.Overlay.read(curvFileName).data
        normed = (curv - np.median(curv))/np.std(curv)

        return normed


def spherePad(surfName, curvFileName, padSize=8, which_norm='Median'):

        surf = fs.Surface.read(surfName)
        curv = normCurvature(curvFileName, which_norm)

        mrisp = surf.parameterize(curv)  # TODO this might be incorrect now
        cols = mrisp.shape[0]
        rows = mrisp.shape[1]

        data = mrisp.squeeze().transpose()

        paddata = np.concatenate((data[rows-padSize:rows, :], data, data[0:padSize, :]), axis=0)

        return paddata
