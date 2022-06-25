import numpy as np
import pyqtgraph as pg
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QApplication, QTreeWidget, QTreeWidgetItem, QHBoxLayout, QWidget


app = QApplication([])
windows = {}


def view_arrays(arrays, *args, **kwargs):
    layers = [{'data': array, 'name': str(i), 'visible': True} for i, array in enumerate(arrays)]
    view(layers, *args, **kwargs)


def view(layers, window_id=0, interactive=False, user_keys_callback=None, title=None,
         legend_width=None, handle_events=True):
    shape = layers[0]['data'].shape
    if not all(layer['data'].shape == shape for layer in layers):
        raise ValueError("All data entries must have the same shape.")
    elif window_id in windows:
        windows[window_id].update_layers(layers)
    elif layers[0]['data'].ndim == 2:
        windows[window_id] = HdavWindow2d(layers, user_keys_callback=user_keys_callback, legend_width=legend_width)
    elif layers[0]['data'].ndim == 3:
        windows[window_id] = HdavWindow3d(layers, user_keys_callback=user_keys_callback, legend_width=legend_width)
    else:
        raise ValueError("All data entries must have 2 or 3 dimensions.")
    if title is not None:
        windows[window_id].setWindowTitle(title)
    windows[window_id].show()
    if handle_events:
        if interactive:
            QApplication.processEvents()
        else:
            QApplication.exec_()

class HdavWindow(QWidget):
    OVERLAY_KEYS = [
        QtCore.Qt.Key_1,
        QtCore.Qt.Key_2,
        QtCore.Qt.Key_3,
        QtCore.Qt.Key_4,
        QtCore.Qt.Key_5,
        QtCore.Qt.Key_6,
        QtCore.Qt.Key_7,
        QtCore.Qt.Key_8,
        QtCore.Qt.Key_9,
        QtCore.Qt.Key_A,
        QtCore.Qt.Key_B,
        QtCore.Qt.Key_C,
        QtCore.Qt.Key_D,
        QtCore.Qt.Key_E,
        QtCore.Qt.Key_F,
        QtCore.Qt.Key_G,
        QtCore.Qt.Key_H,
        QtCore.Qt.Key_I,
        QtCore.Qt.Key_J,
        QtCore.Qt.Key_K,
        QtCore.Qt.Key_L,
        QtCore.Qt.Key_M,
        QtCore.Qt.Key_N,
        QtCore.Qt.Key_O,
        QtCore.Qt.Key_P,
        QtCore.Qt.Key_Q,
        QtCore.Qt.Key_R,
        QtCore.Qt.Key_S,
        QtCore.Qt.Key_T,
        QtCore.Qt.Key_U,
        QtCore.Qt.Key_V,
        QtCore.Qt.Key_W,
        QtCore.Qt.Key_X,
        QtCore.Qt.Key_Y,
        QtCore.Qt.Key_Z,
    ]
    MAX_LAYER_COUNT = len(OVERLAY_KEYS)
    OVERLAY_KEY_LABEL = "123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    def __init__(self, layers, user_keys_callback, legend_width, **kwargs):
        super().__init__(**kwargs)
        if legend_width is None:
            legend_width = 200
        self.layers = layers
        self.user_keys_callback = user_keys_callback or {}
        self.views = pg.GraphicsLayoutWidget(self)
        self.legend = QTreeWidget()
        self.legend.setFixedWidth(legend_width)
        self.legend.setRootIsDecorated(False)
        self.legend.setUniformRowHeights(True)
        self.legend.setAllColumnsShowFocus(True)
        self.legend.setItemsExpandable(False)
        self.legend.setHeaderHidden(True)
        self.legend.setColumnCount(2)
        for i, layer in enumerate(layers[:self.MAX_LAYER_COUNT]):
            item = QTreeWidgetItem([self.OVERLAY_KEY_LABEL[i], layer['name']])
            self.legend.addTopLevelItem(item)
        layout = QHBoxLayout(self)
        layout.addWidget(self.legend)
        layout.addWidget(self.views)

    def keyPressEvent(self, ev):
        key = ev.key()
        if key in self.user_keys_callback:
            callback = self.user_keys_callback[key]
            callback()
            ev.accept()
        elif key in HdavWindow.OVERLAY_KEYS:
            index = HdavWindow.OVERLAY_KEYS.index(key)
            if index < len(self.layers):
                ev.accept()
                self.layers[index]['visible'] = not self.layers[index]['visible']
                self.draw()

    @property
    def shape(self):
        if self.layers:
            return self.layers[0]['data'].shape
        else:
            return None

    def update_layers(self, layers):
        old_shape = self.shape
        old_visibility = [layer['visible'] for layer in self.layers]
        self.layers = layers
        for visibility, layer in zip(old_visibility, self.layers):
            layer['visible'] = visibility
        new_shape = self.shape
        if old_shape != new_shape:
            self.update_cursor(old_shape, new_shape)
        self.draw()

    def draw(self):
        raise NotImplementedError()


class HdavWindow2d(HdavWindow):
    def __init__(self, layers, *args, **kwargs):
        super().__init__(layers, *args, **kwargs)
        self.resize(300, 300)
        self.viewbox = HdavViewBox2d(self)
        self.views.addItem(self.viewbox)

        self.draw()

    def update_cursor(self, old_shape, new_shape):
        pass

    def draw(self):
        self.viewbox.clear()

        for layer in self.layers:
            image = pg.ImageItem(layer['data'])
            if 'cmap' in layer and layer['cmap'] is not None:
                start, stop = layer['cmap'].pos[0], layer['cmap'].pos[-1]
                lut = layer['cmap'].getLookupTable(start=start, stop=stop)
                image.setLookupTable(lut)
            if layer['visible']:
                self.viewbox.addItem(image)


class HdavWindow3d(HdavWindow):
    def __init__(self, layers, *args, **kwargs):
        super().__init__(layers, *args, **kwargs)
        self.resize(900, 300)
        self.viewbox_axial = HdavViewBox3d(self, 0, 1, 2)
        self.viewbox_coronal = HdavViewBox3d(self, 0, 2, 1)
        self.viewbox_sagittal = HdavViewBox3d(self, 1, 2, 0)
        self.views.addItem(self.viewbox_axial)
        self.views.addItem(self.viewbox_coronal)
        self.views.addItem(self.viewbox_sagittal)

        x_pen = pg.mkPen('r')
        y_pen = pg.mkPen('g')
        z_pen = pg.mkPen('b')
        self.cursor_axial_x = pg.InfiniteLine(angle=0, pen=x_pen)
        self.cursor_axial_y = pg.InfiniteLine(angle=90, pen=y_pen)
        self.cursor_coronal_x = pg.InfiniteLine(angle=0, pen=x_pen)
        self.cursor_coronal_z = pg.InfiniteLine(angle=90, pen=z_pen)
        self.cursor_sagittal_y = pg.InfiniteLine(angle=0, pen=y_pen)
        self.cursor_sagittal_z = pg.InfiniteLine(angle=90, pen=z_pen)

        self.cursor_pos = np.array(self.shape) / 2
        self.draw()

    def update_cursor(self, old_shape, new_shape):
        self.set_clipped_cursor(
            [new_shape[dim] * value / old_shape[dim] for dim, value in enumerate(self.cursor_pos)])

    def move_cursor(self, pos):
        self.set_clipped_cursor(pos)
        self.draw()


    def set_clipped_cursor(self, pos):
        min_pos = np.zeros(3, dtype=int)
        max_pos = np.array(self.shape) - 1
        self.cursor_pos = np.maximum(np.minimum(pos, max_pos), min_pos)

    def draw(self):
        self.viewbox_axial.clear()
        self.viewbox_coronal.clear()
        self.viewbox_sagittal.clear()

        for layer in self.layers:
            axial = layer['data'][:, :, int(self.cursor_pos[2])]
            coronal = layer['data'][:, int(self.cursor_pos[1]), :]
            sagittal = layer['data'][int(self.cursor_pos[0]), :, :]
            image_axial = pg.ImageItem(axial)
            image_coronal = pg.ImageItem(coronal)
            image_sagittal = pg.ImageItem(sagittal)
            if 'cmap' in layer and layer['cmap'] is not None:
                start, stop = layer['cmap'].pos[0], layer['cmap'].pos[-1]
                lut = layer['cmap'].getLookupTable(start=start, stop=stop)
                image_axial.setLookupTable(lut)
                image_coronal.setLookupTable(lut)
                image_sagittal.setLookupTable(lut)
            if layer['visible']:
                self.viewbox_axial.addItem(image_axial)
                self.viewbox_coronal.addItem(image_coronal)
                self.viewbox_sagittal.addItem(image_sagittal)

        axial_pos = pg.Point(self.cursor_pos[[0, 1]])
        coronal_pos = pg.Point(self.cursor_pos[[0, 2]])
        sagittal_pos = pg.Point(self.cursor_pos[[1, 2]])
        self.cursor_axial_x.setPos(axial_pos)
        self.cursor_axial_y.setPos(axial_pos)
        self.cursor_coronal_x.setPos(coronal_pos)
        self.cursor_coronal_z.setPos(coronal_pos)
        self.cursor_sagittal_y.setPos(sagittal_pos)
        self.cursor_sagittal_z.setPos(sagittal_pos)
        self.viewbox_axial.addItem(self.cursor_axial_x)
        self.viewbox_axial.addItem(self.cursor_axial_y)
        self.viewbox_coronal.addItem(self.cursor_coronal_x)
        self.viewbox_coronal.addItem(self.cursor_coronal_z)
        self.viewbox_sagittal.addItem(self.cursor_sagittal_y)
        self.viewbox_sagittal.addItem(self.cursor_sagittal_z)


class HdavViewBox(pg.ViewBox):
    def __init__(self, hdav_window, **kwargs):
        kwargs.setdefault('lockAspect', True)
        kwargs.setdefault('border', 'w')
        super().__init__(**kwargs)
        self.hdav_window = hdav_window


class HdavViewBox2d(HdavViewBox):
    def mouseDragEvent(self, ev, axis=None):
        pass

    def wheelEvent(self, ev, axis=None):
        pass


class HdavViewBox3d(HdavViewBox):
    def __init__(self, hdav_window, x_axis, y_axis, z_axis, **kwargs):
        super().__init__(hdav_window, **kwargs)
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.z_axis = z_axis
        self.sensitivity = 0.0002

    def mouseClickEvent(self, ev):
        if ev.button() == QtCore.Qt.LeftButton:
            ev.accept()
            cursor_pos = self.hdav_window.cursor_pos.copy()
            click_pos = self.mapToView(ev.pos())
            cursor_pos[self.x_axis] = click_pos.x()
            cursor_pos[self.y_axis] = click_pos.y()
            self.hdav_window.move_cursor(cursor_pos)

    def mouseDragEvent(self, ev, axis=None):
        if ev.button() == QtCore.Qt.LeftButton:
            ev.accept()
            cursor_pos = self.hdav_window.cursor_pos.copy()
            click_pos = self.mapToView(ev.pos())
            cursor_pos[self.x_axis] = click_pos.x()
            cursor_pos[self.y_axis] = click_pos.y()
            self.hdav_window.move_cursor(cursor_pos)

    def wheelEvent(self, ev, axis=None):
        ev.accept()
        cursor_pos = self.hdav_window.cursor_pos.copy()
        change = ev.delta() * self.sensitivity * self.hdav_window.shape[self.z_axis]
        cursor_pos[self.z_axis] += change
        self.hdav_window.move_cursor(cursor_pos)
