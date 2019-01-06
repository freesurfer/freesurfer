import colorsys
import functools

import numpy as np
from colormath.color_conversions import convert_color
from colormath.color_diff import delta_e_cie2000
from colormath.color_objects import sRGBColor, LabColor

LUMINOSITY_VECTOR = np.array([0.2126, 0.7152, 0.0722])


def hsv_palette(how_many):
    return [colorsys.hsv_to_rgb(segment / how_many, 1.0, 1.0) for segment in range(how_many)]


def luminosity(color):
    return np.dot(LUMINOSITY_VECTOR, np.array(color))


def saturation(color):
    return max(color) - min(color)


def as_gray(r, g, b):
    lum = luminosity([r, g, b])
    return lum, lum, lum


def luminosity_metric(a, b):
    return abs(luminosity(a) - luminosity(b))


def disparity(a, b):
    return disparity_rgb(a[0], a[1], a[2], b[0], b[1], b[2])


@functools.lru_cache(None)
def disparity_rgb(r0, g0, b0, r1, g1, b1):
    color_differencer = cie_2000_differencer
    gray0 = luminosity([r0, g0, b0])
    gray1 = luminosity([r1, g1, b1])
    gray_difference = color_differencer(gray0, gray0, gray0, gray1, gray1, gray1)
    color_difference = color_differencer(r0, g0, b0, r1, g1, b1)
    return color_difference, gray_difference


@functools.lru_cache(None)
def cie_2000_differencer(r0, g0, b0, r1, g1, b1):
    return delta_e_cie2000(lab_color(r0, g0, b0), lab_color(r1, g1, b1))


@functools.lru_cache(None)
def lab_color(r, g, b):
    return convert_color(sRGBColor(r, g, b), LabColor)


def perceptual_difference(a, b):
    color_difference, gray_difference = disparity(a, b)
    # Uses the color vision difference as long as it is not more than twice color blind difference
    return min(color_difference, gray_difference * 2)


def saturation_biased_perceptual_difference(a, b):
    return saturation(a) * perceptual_difference(a, b)


def candidate_color_listing(segments=None):
    if segments is None:
        segments = [10, 40, 4]
    ranges = [[index / segment for index in range(segment + 1)] for segment in segments]
    [red_list, green_list, blue_list] = ranges
    return [[red, green, blue]
            for red in red_list
            for green in green_list
            for blue in blue_list
            ]


def closest_distance(color, comparison_list, metric=None):
    if metric is None:
        metric = saturation_biased_perceptual_difference
    return min([metric(color, other_color) for other_color in comparison_list])


def best_color_index(candidates, comparison_list, metric=None):
    best_color_index = None
    best_distance = -1
    for index, color in enumerate(candidates):
        distance = closest_distance(color, comparison_list, metric)
        if distance > best_distance:
            best_color_index = index
            best_distance = distance
    return best_color_index


def maximal_distance_palette(max_size=35, segments=None, metric=None):
    color_list = candidate_color_listing(segments)
    candidates = color_list[1:-1]
    bad_choices = [color_list[0], color_list[-1]]  # stay away from black and white
    palette = []
    while candidates and len(palette) < max_size:
        comparision_list = bad_choices + palette
        index = best_color_index(candidates, comparision_list, metric)
        best_color = candidates[index]
        palette.append(best_color)
        candidates = candidates[0:index] + candidates[index + 1:]
    return palette


def generate_python_code_for_default_palette(palette):
    print('DEFAULT_PALETTE = [')
    for index, color in enumerate(palette):
        print('    {0},  # luminosity={1}'.format(color, luminosity(color)))
    print(']')


if __name__ == '__main__':
    generate_python_code_for_default_palette(maximal_distance_palette(metric=perceptual_difference))
