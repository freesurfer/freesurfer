import colorsys

import numpy as np
import pyqtgraph as pg

from samseg.hdav import hdav
from samseg.hdav.hdav.hdav import HdavWindow


class _BaseShowFigures:
    def __init__(self, interactive=False, image_alpha=0.6):
        self.interactive = interactive
        self._next_window_id = 0
        self.image_alpha = image_alpha


class DoNotShowFigures(_BaseShowFigures):
    def show(self, **kwargs):
        pass


class ShowFigures(_BaseShowFigures):
    def show(self, images=None, probabilities=None, mesh=None, names=None, window_id=None, shape=None):
        if probabilities is None and mesh:
            if shape is None:
                shape = images.shape
            probabilities = mesh.rasterize_atlas(shape[0:3])
        if probabilities is None:
            probability_layers = []
            image_alpha = None
        else:
            probability_layers = self.probability_layers(probabilities, names)
            image_alpha = self.image_alpha
        if images is None:
            image_layers = []
        else:
            image_layers = self.image_layers(images, alpha=image_alpha)
        probability_max = HdavWindow.MAX_LAYER_COUNT - len(image_layers)
        if len(probability_layers) > probability_max:
            tail = np.sum(layer['data'] for layer in probability_layers[probability_max - 1:])
            probability_layers[probability_max - 1]['data'] = tail
            probability_layers = probability_layers[:probability_max]
        layers = probability_layers + image_layers
        self.show_layers(layers, window_id)

    def probability_layers(self, probabilities, names, prefix=None, visibility=True):
        if prefix is None:
            prefix = 'Label'
        image_list = self.create_image_list(probabilities)
        layer_count = len(image_list)
        if names is None:
            names = self.generate_names(prefix, layer_count)
        palette = hsv_palette(layer_count)
        return [{
            'name': names[index],
            'data': image_list[index],
            'visible': visibility,
            'cmap': transparency_color_map(palette[index])
        } for index in range(layer_count)]

    def image_layers(self, images, prefix=None, visibility=True, alpha=None):
        if prefix is None:
            prefix = 'Contrast'
        image_list = self.create_image_list(images)
        layer_count = len(image_list)
        names = self.generate_names(prefix, layer_count)
        if alpha is None:
            cmap = None
        else:
            start = np.min(images)
            stop = np.max(images)
            cmap = transparency_color_map(
                (1.0, 1.0, 1.0),  # white
                bottom_rgb_color=(0., 0., 0.),
                start=start,
                stop=stop,
                bottom_alpha=alpha,
                top_alpha=alpha)
        return [{
            'name': names[index],
            'data': image_list[index],
            'visible': visibility,
            'cmap': cmap,
        } for index in range(layer_count)]

    def create_image_list(self, source):
        shape = source.shape
        if len(shape) > 3:
            layer_count = shape[3]
            image_list = [source[:, :, :, index] for index in range(layer_count)]
        else:
            image_list = [source]
        return image_list

    def generate_names(self, prefix, count):
        return ['{0} {1}'.format(prefix, index + 1) for index in range(count)]

    def show_layers(self, layers, window_id):
        if window_id is None:
            window_id = str(self._next_window_id)
            self._next_window_id += 1
        if layers:
            hdav.view(layers, interactive=self.interactive, window_id=window_id)


def transparency_color_map(top_rgb_color, bottom_rgb_color=None, start=0, stop=65535, bottom_alpha=0.0, top_alpha=1.0):
    if bottom_rgb_color is None:
        bottom_rgb_color = top_rgb_color
    positions = [start, stop]
    zero_color = list(bottom_rgb_color) + [bottom_alpha]  # fully transparent
    top_color = list(top_rgb_color) + [top_alpha]  # fully opaque
    colors = [zero_color, top_color]
    return pg.ColorMap(positions, colors)


def hsv_palette(how_many):
    return [colorsys.hsv_to_rgb(segment / how_many, 1.0, 1.0) for segment in range(how_many)]
