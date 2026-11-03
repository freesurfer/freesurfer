description = """
This script runs the claustrum segmentation using SynthSeg.
"""

# project imports
from SynthSeg.predict import predict
import numpy as np
import argparse
import sys
import os
import platform

def print_vm_peak():
    """
    Print the VM peak of the running process. This is only available
    on linux platforms.
    """
    if platform.system() != 'Linux':
        return
    procstat = os.path.join('/proc', str(os.getpid()), 'status')
    fp = open(procstat, 'r')
    lines = fp.readlines()
    for line in lines:
        strs = line.split()
        if(len(strs) < 3):
            continue
        if(strs[0] != 'VmPeak:'):
            continue
        print('claustrum-seg VmPeak:', int(strs[1]))


def main():
    default_model_path = os.path.join(os.environ.get('FREESURFER_HOME_FSPYTHON'), 'claustrum_seg_20250616.h5')
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', '--i', help='Image to segment (can be a folder)')
    parser.add_argument('-o', '--o', help='Output seg (can be a folder)');
    parser.add_argument('-c', '--csv', help='Output csv of mm3 volume ');
    parser.add_argument('-p', '--post', help='Output posterior image ');
    parser.add_argument('-m', '--model', default=default_model_path, help='instead of default ');

    # This does nothing
    parser.add_argument("--threads", type=int, default=1, help="(optional) Number of cores to be used. Default is 1.")


    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    # print out the command line
    print(' '.join(sys.argv))
        
    # parse commandline
    args = parser.parse_args()
                                        
    path_images = args.i;
    path_segm = args.o;
    path_posteriors = args.post;
    path_vol = args.csv;
    path_model = args.model


    print('model ',path_model)


    segmentation_labels = np.array([0, 4,15,16,24,41,42,43,44,46,47,49,50,51,52,53,54,58,60,139])

    # We can now provide various parameters to control the preprocessing of the input.
    cropping = None
    target_res = 0.35
    flip = False 
    n_neutral_labels = 20
    sigma_smoothing = 0.5
    topology_classes = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])
    keep_biggest_component = True

    # Regarding the architecture of the network, we must provide the predict function with the same parameters as during
    # training.
    n_levels = 5
    nb_conv_per_level = 2
    conv_size = 3
    unet_feat_count = 24
    activation = 'elu'
    feat_multiplier = 2

    gt_folder = None
    compute_distances = False


    predict(path_images,
            path_segm,
            path_model,
            segmentation_labels,
            n_neutral_labels=n_neutral_labels,
            path_posteriors=path_posteriors,
            path_volumes=path_vol,
            cropping=cropping,
            target_res=target_res,
            flip=flip,
            topology_classes=topology_classes,
            sigma_smoothing=sigma_smoothing,
            keep_biggest_component=keep_biggest_component,
            n_levels=n_levels,
            nb_conv_per_level=nb_conv_per_level,
            conv_size=conv_size,
            unet_feat_count=unet_feat_count,
            feat_multiplier=feat_multiplier,
            activation=activation,
            gt_folder=gt_folder,
            compute_distances=compute_distances)

    print_vm_peak()
    print('Claustrum seg done')


# execute script
if __name__ == '__main__':
    main()
