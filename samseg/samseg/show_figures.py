import math

import numpy as np
import pyqtgraph as pg
from PyQt5 import QtGui

from samseg.hdav import hdav
from samseg.hdav.hdav.hdav import HdavWindow

DEFAULT_PALETTE = [
    [0.7, 0.15, 1.0],  # luminosity=0.3283
    [0.3, 0.65, 0.0],  # luminosity=0.5286599999999999
    [1.0, 0.725, 0.0],  # luminosity=0.73112
    [0.7, 0.05, 0.0],  # luminosity=0.18458000000000002
    [0.6, 0.375, 0.5],  # luminosity=0.43186
    [0.4, 0.65, 1.0],  # luminosity=0.62212
    [0.4, 0.975, 1.0],  # luminosity=0.85456
    [0.2, 0.2, 1.0],  # luminosity=0.25776
    [0.0, 0.1, 0.5],  # luminosity=0.10762
    [0.7, 0.3, 0.25],  # luminosity=0.38143
    [0.4, 0.45, 1.0],  # luminosity=0.47907999999999995
    [0.9, 0.6, 0.75],  # luminosity=0.67461
    [0.9, 1.0, 0.25],  # luminosity=0.9245899999999999
    [0.1, 1.0, 0.75],  # luminosity=0.7906099999999999
    [0.1, 0.175, 0.0],  # luminosity=0.14642
    [0.2, 0.275, 0.75],  # luminosity=0.29335
    [0.2, 0.225, 0.25],  # luminosity=0.22149
    [0.0, 0.7, 1.0],  # luminosity=0.57284
    [0.3, 0.0, 0.0],  # luminosity=0.06378
    [0.4, 0.425, 0.25],  # luminosity=0.40704999999999997
    [1.0, 0.2, 0.0],  # luminosity=0.35564
    [0.6, 0.475, 0.5],  # luminosity=0.5033799999999999
    [0.8, 0.4, 0.0],  # luminosity=0.45616
    [0.2, 0.75, 0.25],  # luminosity=0.59697
    [0.1, 0.775, 1.0],  # luminosity=0.64774
    [0.9, 0.925, 0.5],  # luminosity=0.889
    [0.4, 0.6, 0.5],  # luminosity=0.55026
    [0.8, 0.8, 0.25],  # luminosity=0.76029
    [0.5, 0.925, 0.75],  # luminosity=0.82201
    [0.9, 0.975, 1.0],  # luminosity=0.96086
    [0.7, 0.775, 0.0],  # luminosity=0.7031000000000001
    [0.2, 0.35, 0.25],  # luminosity=0.31089
    [0.2, 0.2, 0.75],  # luminosity=0.23971
    [0.6, 0.0, 0.0],  # luminosity=0.12756
    [0.7, 0.0, 0.75],  # luminosity=0.20297
]


class DoNotShowFigures:
    def show(self, **kwargs):
        pass

    def plot(self, **kwargs):
        pass


class ShowFigures:
    def __init__(self, interactive=False, image_alpha=0.6, palette=None):
        if palette is None:
            palette = DEFAULT_PALETTE
        self.palette = palette
        self.interactive = interactive
        self._next_window_id = 0
        self.image_alpha = image_alpha

    def create_palette(self, layer_count):
        palette = []
        self_size = len(self.palette)
        while layer_count > self_size:
            palette += self.palette
            layer_count -= self_size
        if layer_count:
            palette += self.palette[:layer_count]
        return palette

    def show(self, auto_scale=False, images=None, image_list=None, probabilities=None, mesh=None, names=None,
             window_id=None, shape=None, title=None):
        if image_list is None:
            image_list = []
        if images is not None:
            image_list += self.create_image_list(images)
        if probabilities is None and mesh:
            if shape is None:
                shape = images.shape
            probabilities = mesh.rasterize_atlas(shape[0:3])
        if probabilities is None:
            probability_layers = []
            image_alpha = 1.0 if auto_scale else None
        else:
            probability_layers = self.probability_layers(probabilities, names)
            image_alpha = self.image_alpha
        image_layers = self.image_layers(image_list, alpha=image_alpha)
        probability_max = HdavWindow.MAX_LAYER_COUNT - len(image_layers)
        if len(probability_layers) > probability_max:
            tail = np.sum(layer['data'] for layer in probability_layers[probability_max - 1:])
            probability_layers[probability_max - 1]['data'] = tail
            probability_layers = probability_layers[:probability_max]
        layers = probability_layers + image_layers
        self.show_layers(layers, window_id, title)

    def probability_layers(self, probabilities, names, prefix=None, visibility=True):
        if prefix is None:
            prefix = 'Label'
        image_list = self.create_image_list(probabilities)
        layer_count = len(image_list)
        if names is None:
            names = self.generate_names(prefix, layer_count)
        palette = self.create_palette(layer_count)
        return [{
            'name': names[index],
            'data': image_list[index],
            'visible': visibility,
            'cmap': transparency_color_map(palette[index])
        } for index in range(layer_count)]

    def image_layers(self, image_list, prefix=None, visibility=True, alpha=None):
        if len(image_list) == 0:
            return []
        if prefix is None:
            prefix = 'Contrast'
        layer_count = len(image_list)
        names = self.generate_names(prefix, layer_count)
        if alpha is None:
            cmap = None
        else:
            start = min([np.min(image) for image in image_list])
            stop = max([np.max(image) for image in image_list])
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
        if source is None:
            return []
        shape = source.shape
        if len(shape) > 3:
            layer_count = shape[3]
            image_list = [source[:, :, :, index] for index in range(layer_count)]
        else:
            image_list = [source]
        return image_list

    def generate_names(self, prefix, count):
        return ['{0} {1}'.format(prefix, index + 1) for index in range(count)]

    def show_layers(self, layers, window_id, title):
        if window_id is None:
            window_id = str(self._next_window_id)
            self._next_window_id += 1
        if layers:
            hdav.view(layers, interactive=self.interactive, window_id=window_id, title=title)

    def plot(self, data, title=None):
        pg.plot(y=data, pen=pg.mkPen(color='l', width=3), title=title)
        if self.interactive:
            QtGui.QApplication.processEvents()
        else:
            QtGui.QApplication.exec_()


def transparency_color_map(top_rgb_color, bottom_rgb_color=None, start=0, stop=65535, bottom_alpha=0.0, top_alpha=1.0):
    if bottom_rgb_color is None:
        bottom_rgb_color = top_rgb_color
    positions = [start, stop]
    zero_color = list(bottom_rgb_color) + [bottom_alpha]  # fully transparent
    top_color = list(top_rgb_color) + [top_alpha]  # fully opaque
    colors = [zero_color, top_color]
    return pg.ColorMap(positions, colors)


def show_palette(palette=None):
    if palette is None:
        palette = DEFAULT_PALETTE
    OFF = 0
    ON = 65535
    color_count = len(palette)
    unit = 5
    stripe = 3 * unit
    width = stripe * color_count
    image = np.zeros([width, width, unit, color_count])

    # make horizontal stripes
    for color in range(color_count):
        base_row = color * stripe
        end_row = base_row + stripe
        image[:, base_row:end_row, :, color] = ON

        # with inserted squares
        for square_color in range(color_count):
            square_base_row = base_row + unit
            square_end_row = square_base_row + unit
            square_base_column = square_color * stripe + unit
            square_end_column = square_base_column + unit
            image[square_base_column:square_end_column, square_base_row:square_end_row, :, color] = OFF
            image[square_base_column:square_end_column, square_base_row:square_end_row, :, square_color] = ON
    ShowFigures(palette=palette).show(probabilities=image)


def show_example_plot():
    data = [math.sqrt(val) for val in range(10)]
    ShowFigures().plot(data, title='Hello')


if __name__ == '__main__':
    # show_palette()
    show_example_plot()
