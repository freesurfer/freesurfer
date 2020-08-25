import os
import sys
from argparse import ArgumentParser
import logging
import tensorflow as tf
tf.get_logger().setLevel(logging.ERROR)

hypothalamus_home = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))
sys.path.append(hypothalamus_home)
from hypothalamus_seg.predict import predict

parser = ArgumentParser()

# positional arguments
parser.add_argument("path_images", type=str, help="images to segment. Can be the path to a single image or to a folder")
parser.add_argument("path_segmentations", type=str, help="path where to save the segmentations. Must be the same type "
                                                         "as path_images (path to a single image or to a folder)")

# preprocessing/postprocessing parameters
parser.add_argument("--out_posteriors", type=str, default=None, dest="path_posteriors",
                    help="path where to save the posteriors. Must be the same type as path_images (path to a single "
                         "image or to a folder)")
parser.add_argument("--out_volumes", type=str, default=None, dest="path_volumes",
                    help="path to a csv file where to save the volumes of all subunits for all patients")

args = vars(parser.parse_args())
args['path_model'] = os.path.join(hypothalamus_home, 'data', 'model.h5')
predict(**args)
