from samseg.show_figures import ShowFigures


class TestShowFigures:
    def setup(self):
        self.palette = [
            [1,2,3],
            [4,5,6],
            [7,8,9],
        ]
        self.show_figures = ShowFigures(palette=self.palette)

    def test_create_palette(self):
        for layer_count in [2, 3, 5]:
            palette = self.show_figures.create_palette(layer_count)
            assert len(palette) == layer_count
            for index in range(layer_count):
                assert palette[index] == self.palette[index % 3]