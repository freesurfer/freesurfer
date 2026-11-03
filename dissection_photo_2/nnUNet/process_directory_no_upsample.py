import torch
from argparse import ArgumentParser
import os
import glob
import shutil
import cv2
import sys
sys.path.insert(0, os.path.split(__file__)[0])
from nnunetv2.inference.predict_from_raw_data import nnUNetPredictor

def main(args):

    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    if os.path.isdir(args.out_dir) is False:
        os.mkdir(args.out_dir)
    tempdir = args.out_dir + '/tmp/'
    if os.path.isdir(tempdir):
        shutil.rmtree(tempdir, ignore_errors=True)
    os.mkdir(tempdir)

    # Same formats as the GUI's ImageNameFilters(), in any case
    g = sorted(f for f in glob.glob(os.path.join(glob.escape(args.in_dir), '*'))
               if os.path.splitext(f)[1].lower() in ('.jpg', '.jpeg', '.png', '.tif', '.tiff', '.bmp'))

    print('Downsampling ' + str(len(g)) + ' files')
    for input in g:
        I = cv2.imread(input)
        Id = cv2.resize(I, [0,0], fx = 0.2, fy = 0.2, interpolation=cv2.INTER_AREA)
        # Name up to the last dot, like the GUI's completeBaseName(), so .jpeg/.tiff work too
        name = os.path.splitext(os.path.basename(input))[0]
        cv2.imwrite(os.path.join(tempdir, name + '_0000.png'), Id[:, :, 2])
        cv2.imwrite(os.path.join(tempdir, name + '_0001.png'), Id[:, :, 1])
        cv2.imwrite(os.path.join(tempdir, name + '_0002.png'), Id[:, :, 0])

    print('Loading model')
    perform_everything_on_device = True if device=='cuda' else False
    predictor = nnUNetPredictor(device=torch.device(device), \
                                perform_everything_on_device=perform_everything_on_device, \
                                verbose=False, verbose_preprocessing=False)

    # One fold on the CPU: speed, and multi-fold prediction crashes there
    if device=='cpu':
        predictor.initialize_from_trained_model_folder(args.model_path, use_folds=(0,))
    else:
        predictor.initialize_from_trained_model_folder(args.model_path, use_folds=(0, 1, 2, 3, 4))

    print('Predicting!')
    predictor.predict_from_files(tempdir, args.out_dir, save_probabilities=True)

    print('Cleaning up')
    shutil.rmtree(tempdir, ignore_errors=True)

    print('All done')


if __name__ == "__main__":
    parser = ArgumentParser()

    parser.add_argument("--in_dir", type=str, dest="in_dir", default=None)
    parser.add_argument("--out_dir", type=str, dest="out_dir", default=None)
    parser.add_argument("--model_path", type=str, dest="model_path", default=None)

    args = parser.parse_args()

    main(args)
