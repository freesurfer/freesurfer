#! /usr/bin/env python


import argparse
import os

# from psacnn_brain_segmentation import predict_segmentation

parser = argparse.ArgumentParser(
    description='''sscnn_skullstrip
            Skull strip an infant MRI volume
            -i or --input-file: Input volume nii or mgz file path
            -o or --output-dir: A filepath to the output directory. 
            -c or --contrast('t1w' supported)
            -gpu: GPU id to be used for computation

                ''')
parser.add_argument('-i', '--input-file', type=str, required=True)
parser.add_argument('-o', '--output-dir', type=str, required=True)
parser.add_argument('-c', '--contrast', type=str, required=False)
parser.add_argument('-gpu', nargs='?', const=0, default=0, type=int)
args = parser.parse_args()
print(args)
for argname in ['input_file']:
    setattr(args, argname, os.path.abspath(os.path.expanduser(getattr(args, argname))))
    if not os.path.exists(getattr(args, argname)):
        raise RuntimeError('%s: %s does not exist. Exiting.' % (argname.replace('_', ' ').title(),
                                                                getattr(args, argname)))



from preprocess import sscnn_workflow

print(args.gpu)

if args.gpu < 0 :
    use_gpu = False
    batch_size = 16
    sample_rate = 40000
else:
    use_gpu = True
    batch_size=4
    sample_rate = 20000


# when run as a workflow save_label_image = False (as the nipype workflow will save it
    # input_dir =  '/autofs/space/vault_007/users/lzollei/NadineGaab/Work/new_skullstrips-hf/poor_skullstrip/INF045.hf/'
                 # /autofs/space/vault_007/users/lzollei/AmodJog/data/NadineGaab-Bangladesh'
    # output_dir = '/autofs/space/vault_007/users/lzollei/AmodJog/code/debug/'
        # '/autofs/space/vault_007/users/lzollei/AmodJog/results/NadineGaab-Bangladesh/'
    # gpu_id = 0
sscnn_workflow(input_file=args.input_file,
                output_dir=args.output_dir,
                contrast=args.contrast,
                use_gpu=use_gpu,
                gpu_id=args.gpu,
                save_label_image=False,
                save_prob_image=False,
                batch_size=batch_size,
               )


# from psacnn_brain_segmentation import predict
# predict.predict_segmentation(input_file=args.input_file, output_dir=args.output_dir, contrast=args.contrast,
#                              model_file=args.model_file, output_soft=args.output_soft,patch_dim=args.patch_dim)

