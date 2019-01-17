from .bindings import CoreSurface


class Surface(CoreSurface):
    def __init__(self, filename):
        super(Surface, self).__init__(filename)
