class ColorScheme:
    def __init__(self, number_of_classes, gm_class_number):
        self.gm_class_number = gm_class_number
        self.number_of_classes = number_of_classes
        self.colors = generate_colors(number_of_classes)

def generate_colors(number_of_classes):
    # colors = 255 * [ hsv( numberOfClasses ) ones( numberOfClasses, 1 ) ];
    pass