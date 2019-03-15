#! /usr/bin/env python


import argparse
import os

# from psacnn_brain_segmentation import predict_segmentation

parser = argparse.ArgumentParser(
    description='''psacnn_segment
            Segment a skullstripped, inhomogeneity corrected, brain MRI images in the FreeSurfer conformed space.
            -i or --input-file: Input MRI (skullstripped, inhomogeneity-corrected) nii or mgz file path
            -o or --output-dir: A filepath to the output directory. 
            -m or --model-file: model file to be used
            -gpu: GPU id to be used for computation
            --output-soft: Probabilistic output image is saved
            -p or --patch-dim: patch size to be used

                ''')
parser.add_argument('-i', '--input-file', type=str, required=True)
parser.add_argument('-o', '--output-dir', type=str, required=True)
parser.add_argument('-pre', '--preprocess', action='store_true')
parser.add_argument('-m', '--model-file', type=str, required=False)
parser.add_argument('-c', '--contrast', type=str, required=False)
parser.add_argument('-gpu', nargs='?', const=0, default=0, type=int)
# parser.add_argument('-os', '--output-soft', action='store_true')
parser.add_argument('-p', '--patch-size', type=int, choices=[64, 96], default=96)
# parser.add_argument('-b', '--batch-size', type=int, required=False, default=4)
args = parser.parse_args()
print(args)
for argname in ['input_file']:
    setattr(args, argname, os.path.abspath(os.path.expanduser(getattr(args, argname))))
    if not os.path.exists(getattr(args, argname)):
        raise RuntimeError('%s: %s does not exist. Exiting.' % (argname.replace('_', ' ').title(),
                                                                getattr(args, argname)))





from preprocess import psacnn_workflow

#print(args.gpu)

if args.gpu < 0 :
    use_gpu = False
    batch_size = 16
    sample_rate = 40000
else:
    use_gpu = True
    batch_size=4
    sample_rate = 20000


# when run as a workflow save_label_image = False (as the nipype workflow will save it
psacnn_workflow(input_file=args.input_file,
                output_dir=args.output_dir,
                use_preprocess=args.preprocess,
                model_file=args.model_file,
                contrast=args.contrast,
                use_gpu=use_gpu,
                gpu_id=args.gpu,
                save_label_image=False,
                save_prob_image=False,
                patch_size=args.patch_size,
                batch_size=batch_size,
                sample_rate=sample_rate)

# from psacnn_brain_segmentation import predict
# predict.predict_segmentation(input_file=args.input_file, output_dir=args.output_dir, contrast=args.contrast,
#                              model_file=args.model_file, output_soft=args.output_soft,patch_dim=args.patch_dim)

