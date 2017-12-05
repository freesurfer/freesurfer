from optparse import OptionParser

from easydict import EasyDict


def parse_args(argv=None, parser=OptionParser()):
    parser.add_option("-o", "--output", metavar="FOLDER", help="output to FOLDER")
    parser.add_option('-i', '--input', action='append', metavar="FILE", dest='image_file_names',
                      help="input image(s) from FILE")
    parser.add_option('--threads', type='int', default=1, help="number of threads")
    parser.add_option("-r", "--regmat", metavar="FILE", help="skip registration and read from FILE")
    parser.add_option('-m', '--missing', dest='missing_structures', action='append', metavar="LABEL",
                      help="LABEL is a missing structure (repeat for multiple missing labels)")
    parser.add_option('--exvivo', action='store_true', default=False, help="run samseg exvivo")
    parser.add_option('-v', '--verbose', action='store_true', default=False, help="verbose debug output")
    options, args = parser.parse_args(args=argv)
    if options.missing_structures is None:
        options.missing_structures = []
    if not options.image_file_names:
        parser.error("must specify at least one input")
    return EasyDict(options.__dict__)

if __name__ == '__main__':
    print(parse_args())