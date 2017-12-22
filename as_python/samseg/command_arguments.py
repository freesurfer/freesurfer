import argparse


def parse_args(argv=None, parser=argparse.ArgumentParser()):
    parser.add_argument("-o", "--output", metavar="FOLDER", help="output to FOLDER")
    parser.add_argument('-i', '--input', action='append', metavar="FILE", dest='image_file_names',
                      help="input image(s) from FILE")
    parser.add_argument('--threads', type=int, default=1, help="number of threads")
    parser.add_argument("-r", "--regmat", metavar="FILE", help="skip registration and read from FILE")
    parser.add_argument('-m', '--missing', dest='missing_structures', action='append', metavar="LABEL",
                      help="LABEL is a missing structure (repeat for multiple missing labels)")
    parser.add_argument('--showfigs', action='store_true', default=False, help="show figures during run")
    parser.add_argument('--nobrainmask', action='store_true', default=False, help="no initial brain masking based on affine atlas registration")
    parser.add_argument('--diagcovs', action='store_true', default=False, help="use diagonal covariance matrices (only affect multi-contrast case)")
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help="verbose debug output")
    args = parser.parse_args(args=argv)
    if not args.image_file_names:
        parser.error("must specify at least one input")
    return args

if __name__ == '__main__':
    print(parse_args())