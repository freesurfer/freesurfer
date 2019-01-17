from .bindings import CoreVolume


class Volume(CoreVolume):
    def __init__(self, filename):
        super(Volume, self).__init__(filename)