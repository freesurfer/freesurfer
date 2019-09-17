import numpy as np


class Colormap:

    def __init__(self, colorlist):
        self._segments = np.linspace(0, 1, len(colorlist))
        self._colors = np.array(colorlist)

    def __call__(self, samples):
        samples = np.array(samples)
        indices = np.argmax(samples[:, np.newaxis] < self._segments, axis=1) - 1
        indices[indices < 0] = len(self._segments) - 2
        return self._retrieve(samples, indices).astype('int')

    def list(self, num):
        raise NotImplementedError()

    def _retrieve(self, samples, indices):
        raise NotImplementedError()


class ListedColormap(Colormap):

    def _retrieve(self, samples, indices):
        return self._colors[indices]

    def list(self, num):
        return np.resize(self._colors, (num, 3))


class LinearColormap(Colormap):

    def _retrieve(self, samples, indices):
        low = self._segments[indices]
        high = self._segments[indices + 1]
        ratio = ((samples - low) / (high - low))[:, np.newaxis]
        return self._colors[indices + 1] * ratio + self._colors[indices] * (1 - ratio)

    def list(self, num):
        return self(np.linspace(0, 1, num))


cmaps = {}


def get_cmap(cmap):
    if isinstance(cmap, str):
        return cmaps.get(cmap)
    elif isinstance(cmap, Colormap):
        return cmap
    else:
        return None


cmaps['spectral'] = LinearColormap([
    [158,   1,  66],
    [213,  62,  79],
    [244, 109,  67],
    [253, 174,  97],
    [254, 224, 139],
    [255, 255, 191],
    [230, 245, 152],
    [171, 221, 164],
    [102, 194, 165],
    [ 50, 136, 189],
    [ 94,  79, 162]
])

cmaps['pastel'] = ListedColormap([
    [141, 211, 199],
    [255, 255, 179],
    [190, 186, 218],
    [251, 128, 114],
    [128, 177, 211],
    [253, 180,  98],
    [179, 222, 105],
    [252, 205, 229],
    [217, 217, 217],
    [188, 128, 189],
    [204, 235, 197],
    [255, 237, 111]
])
