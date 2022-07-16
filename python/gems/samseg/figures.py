import math
import numpy as np


# Color palette generated via main method of ColorScheme.py
DEFAULT_PALETTE = [
    [0.7, 0.15, 1.0],    # luminosity=0.32830
    [0.3, 0.65, 0.0],    # luminosity=0.52866
    [1.0, 0.725, 0.0],   # luminosity=0.73112
    [0.7, 0.05, 0.0],    # luminosity=0.18458
    [0.6, 0.375, 0.5],   # luminosity=0.43186
    [0.4, 0.65, 1.0],    # luminosity=0.62212
    [0.4, 0.975, 1.0],   # luminosity=0.85456
    [0.2, 0.2, 1.0],     # luminosity=0.25776
    [0.0, 0.1, 0.5],     # luminosity=0.10762
    [0.7, 0.3, 0.25],    # luminosity=0.38143
    [0.4, 0.45, 1.0],    # luminosity=0.47908
    [0.9, 0.6, 0.75],    # luminosity=0.67461
    [0.9, 1.0, 0.25],    # luminosity=0.92459
    [0.1, 1.0, 0.75],    # luminosity=0.79061
    [0.1, 0.175, 0.0],   # luminosity=0.14642
    [0.2, 0.275, 0.75],  # luminosity=0.29335
    [0.2, 0.225, 0.25],  # luminosity=0.22149
    [0.0, 0.7, 1.0],     # luminosity=0.57284
    [0.3, 0.0, 0.0],     # luminosity=0.06378
    [0.4, 0.425, 0.25],  # luminosity=0.40705
    [1.0, 0.2, 0.0],     # luminosity=0.35564
    [0.6, 0.475, 0.5],   # luminosity=0.50338
    [0.8, 0.4, 0.0],     # luminosity=0.45616
    [0.2, 0.75, 0.25],   # luminosity=0.59697
    [0.1, 0.775, 1.0],   # luminosity=0.64774
    [0.9, 0.925, 0.5],   # luminosity=0.88900
    [0.4, 0.6, 0.5],     # luminosity=0.55026
    [0.8, 0.8, 0.25],    # luminosity=0.76029
    [0.5, 0.925, 0.75],  # luminosity=0.82201
    [0.9, 0.975, 1.0],   # luminosity=0.96086
    [0.7, 0.775, 0.0],   # luminosity=0.70310
    [0.2, 0.35, 0.25],   # luminosity=0.31089
    [0.2, 0.2, 0.75],    # luminosity=0.23971
    [0.6, 0.0, 0.0],     # luminosity=0.12756
    [0.7, 0.0, 0.75],    # luminosity=0.20297
]
DEFAULT_PALETTE = ( np.array( DEFAULT_PALETTE ) * 255 ).round().astype( 'int' ).tolist()



def import_graphical_libraries():
    try:
        global pg, QApplication, view, HdavWindow
        import pyqtgraph as pg
        from PyQt5.QtWidgets import QApplication
        
        from .hdav import view, HdavWindow
    except ImportError:
       return False
    else:
       return True


def initVisualizer(showfigs, movie):
    if showfigs or movie:
        if not import_graphical_libraries():
            raise RuntimeError('the samseg visualization tool requires pyqtgraph and pyqt5 - '
                'please install these python packages to show figures')
        else:
            return ShowFigures(show_flag=showfigs, movie_flag=movie, interactive=True)
    else:
        return DoNotShowFigures()


class DoNotShowFigures:
    def __init__(self, **kwargs):
        pass

    def __str__(self):
        return 'No visualization'

    def show(self, **kwargs):
        pass

    def plot(self, what, **kwargs):
        pass

    def start_movie(self, **kwargs):
        pass

    def show_movie(self, **kwargs):
        pass


class ShowFigures:
    def __str__(self):
        visualization = 'Show'
        if self.movie_flag:
            visualization += ' movies'
        if self.show_flag:
            visualization += ' plots, and images'
        return visualization

    def __init__(self, interactive=False, image_alpha=153, palette=None, show_flag=False, movie_flag=False):
        if palette is None:
            palette = DEFAULT_PALETTE
        self.palette = palette
        self.interactive = interactive
        self._next_window_id = 0
        self.image_alpha = image_alpha
        self.movie_flag = movie_flag
        self.show_flag = show_flag
        self.movies = {}
        self.user_callbacks = {}
        self.plotWindows = {}
        self.nextPlotWindow_id = 0


    def create_palette(self, layer_count):
        self_size = len(self.palette)
        return self.palette * (layer_count // self_size) + self.palette[:layer_count % self_size]

    def show(self, auto_scale=False, images=None, image_list=None, probabilities=None, mesh=None, names=None,
             window_id=None, shape=None, title=None, legend_width=None):
        if image_list is None:
            image_list = []
        if images is not None:
            image_list += self.create_image_list(images)
        if probabilities is None and mesh is not None:
            if shape is None:
                shape = images.shape
            probabilities = mesh.rasterize_atlas(shape[0:3])
        if probabilities is None:
            probability_layers = []
            image_alpha = 255 if auto_scale else None
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
        self.save_movie_frame(window_id, layers)
        if self.show_flag:
            self.show_layers(layers, window_id, title, legend_width)

    def probability_layers(self, probabilities, names, prefix=None, visibility=True):
        if prefix is None:
            prefix = 'Label'
        image_list = self.create_image_list(probabilities)
        layer_count = len(image_list)
        if names is None:
            names = self.generate_names(prefix, layer_count)
        palette = self.create_palette(layer_count)
        return [{
            'name': name,
            'data': image,
            'visible': visibility,
            'cmap': transparency_color_map(layer_color),
        } for name, image, layer_color in zip(names, image_list, palette)]

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
                (255, 255, 255),  # white
                bottom_rgb_color=(0, 0, 0),
                start=start,
                stop=stop,
                bottom_alpha=alpha,
                top_alpha=alpha)
        return [{
            'name': name,
            'data': image,
            'visible': visibility,
            'cmap': cmap,
        } for name, image in zip(names, image_list)]

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

    def show_layers(self, layers, window_id, title, legend_width=None):
        if window_id is None:
            window_id = str(self._next_window_id)
            self._next_window_id += 1
        self.hdav_view(layers, window_id, title, legend_width=legend_width)

    def hdav_view(self, layers, window_id, title, handle_events=True, legend_width=None):
        if layers:
            view(layers,
                 interactive=self.interactive,
                 window_id=window_id,
                 title=title,
                 user_keys_callback=self.user_callbacks.get(window_id),
                 handle_events=handle_events,
                 legend_width=legend_width,
            )

    def save_movie_frame(self, window_id, layers):
        if self.movie_flag:
            layer_sequence = self.movies.get(window_id)
            if layer_sequence:
                # This point can only be reached if a movie was started with this window_id
                layer_sequence.add(layers)

    def start_movie(self, window_id, title=None):
        if self.movie_flag:
            layer_sequence = LayerSequence(title)
            self.movies[window_id] = layer_sequence

            def update_movie(layers):
                # Call backs should not handle events as they are inside an event loop already
                self.hdav_view(layers, window_id, layer_sequence.title, handle_events=False)

            def backward():
                update_movie(layer_sequence.previous())

            def forward():
                update_movie(layer_sequence.next())

            def beginning_of_movie():
                update_movie(layer_sequence.rewind())

            def end_of_movie():
                update_movie(layer_sequence.skip_to_end())

            self.user_callbacks[window_id] = {
                pg.QtCore.Qt.Key_Left: backward,
                pg.QtCore.Qt.Key_Right: forward,
                pg.QtCore.Qt.Key_Up: beginning_of_movie,
                pg.QtCore.Qt.Key_Down: end_of_movie,
            }

    def show_movie(self, window_id, clear_memory=True):
        if self.movie_flag:
            layer_sequence = self.movies.get(window_id)
            if layer_sequence:
                #layers = layer_sequence.rewind()
                layers = layer_sequence.skip_to_end()
                title = layer_sequence.title
                self.hdav_view(layers, window_id, title)
                if clear_memory:
                    self.movies[window_id] = None
                    self.user_callbacks[window_id] = None

    def plot(self, data, title=None, window_id=None ):
        if self.show_flag:
            if window_id is None:
                window_id = str( self.nextPlotWindow_id )
                self.nextPlotWindow_id += 1
            if window_id not in self.plotWindows:
                self.plotWindows[ window_id ] = pg.PlotWidget()
            self.plotWindows[ window_id ].clear()
            self.plotWindows[ window_id ].plot( y=data, pen=pg.mkPen( color='l', width=3 ) )
            if title is not None:
                self.plotWindows[ window_id ].setWindowTitle( title )
            self.plotWindows[ window_id ].show()
            if self.interactive:
                QApplication.processEvents()
            else:
                QApplication.exec_()

class LayerSequence:
    def __init__(self, title):
        self._title = title
        self.time_index = -1
        self.history = []

    @property
    def title(self):
        return '{0} {1}/{2}'.format(self._title, self.time_index + 1, self.frame_count)

    @property
    def current_frame(self):
        if self.time_index >= 0:
            return self.history[self.time_index]
        else:
            return None

    @property
    def frame_count(self):
        return len(self.history)

    def add(self, layers):
        for layer in layers:
            layer['data'] = np.copy(layer['data'])
        self.history.append(layers)
        self.time_index = self.frame_count - 1

    def rewind(self):
        if self.history:
            self.time_index = 0
        else:
            self.time_index = -1
        return self.current_frame

    def skip_to_end(self):
        if self.history:
            self.time_index = self.frame_count - 1
        else:
            self.time_index = -1
        return self.current_frame

    def next(self):
        if self.time_index + 1 < self.frame_count:
            self.time_index += 1
        return self.current_frame

    def previous(self):
        if self.time_index > 0:
            self.time_index -= 1
        return self.current_frame


def transparency_color_map(top_rgb_color, bottom_rgb_color=None, start=0, stop=65535, bottom_alpha=0, top_alpha=255):
    if bottom_rgb_color is None:
        bottom_rgb_color = top_rgb_color
    positions = [start, stop]
    zero_color = list(bottom_rgb_color) + [bottom_alpha]  # fully transparent
    top_color = list(top_rgb_color) + [top_alpha]  # fully opaque
    colors = [zero_color, top_color]
    return pg.ColorMap(positions, colors)


if __name__ == '__main__':
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
        ShowFigures(palette=palette, show_flag=True).show(probabilities=image)


    def show_example_plot():
        data = [math.sqrt(val) for val in range(10)]
        ShowFigures(show_flag=True).plot(data, title='Hello')


    show_palette()
    show_example_plot()
