import argparse

def main():

    parser = argparse.ArgumentParser(description="WMH-SynthSeg: joint segmentation of anatomy and white matter hyperintensities ", epilog='\n')
    parser.add_argument("--i", help="Input image or directory.", required=True)
    parser.add_argument("--o", help="Output segmentation (or directory, if the input is a directory)", required=True)
    parser.add_argument("--csv_vols", help="(optional) CSV file with volumes of ROIs")
    parser.add_argument("--device", default='cpu', help="device (cpu or cuda; optional)")
    parser.add_argument("--threads", type=int, default=1,
                        help="(optional) Number of CPU cores to be used. Default is 1. You can use -1 to use all available cores")

    args = parser.parse_args()

    input_path = args.i
    output_path = args.o
    output_csv_path = args.csv_vols
    device = args.device
    threads = args.threads

    # Prepare list of images to segment and leave before loading packages if nothing to do
    import os
    import sys
    if os.path.exists(input_path) is False:
        raise Exception('Input does not exist')

    if output_csv_path is not None:
        head, tail = os.path.split(output_csv_path)
        if ((len(head)>0) and (os.path.isdir(head) is False)):
            raise Exception('Parent directory of CSV file does not exist')
        if tail.endswith('.csv') is False:
            raise Exception('CSV output must be a CSV file')

    if os.path.isfile(input_path): # file
        if (input_path.endswith('.nii') or input_path.endswith('.nii.gz') or input_path.endswith('.mgz')) is False:
            raise Exception('Input image is not of a supported type (.nii, .nii.gz, or .mgz)')
        head, tail = os.path.split(output_path)
        if len(tail) == 0:
            raise Exception('If input is a file, output must be a file')
        if (tail.endswith('.nii') or tail.endswith('.nii.gz') or tail.endswith('.mgz')) is False:
            raise Exception('Output image is not of a supported type (.nii, .nii.gz, or .mgz)')
        if ((len(head) > 0) and (os.path.isdir(head) is False)):
            raise Exception('Parent directory of output image does not exist')
        images_to_segment = [input_path]
        segmentations_to_write = [output_path]

    if os.path.isdir(input_path):  # directory
        images_to_segment = []
        segmentations_to_write = []
        for im in os.listdir(input_path):
            if im.endswith('.nii') or im.endswith('.nii.gz') or im.endswith('.mgz'):
                images_to_segment.append(os.path.join(input_path, im))
                segmentations_to_write.append(os.path.join(output_path, im.replace('.nii', '_seg.nii').replace('.mgz', '_seg.mgz')))
        if len(images_to_segment) == 0:
            raise Exception('Input directory does not contain images with supported type (.nii, .nii.gz, or .mgz)')
        if output_path.endswith('.nii') or output_path.endswith('.nii.gz') or output_path.endswith('.mgz'):
            raise Exception('If input is a directory, output should be a directory too')
        if os.path.isdir(output_path) is False:
            os.mkdir(output_path)

    # We only import packages if we managed to parse
    print('Arguments seem correct; loading Python packages...')
    import torch
    sys.path.append(os.path.abspath(os.path.dirname(__file__)))
    from unet3d.model import UNet3D
    from utils import MRIread, MRIwrite, myzoom_torch, align_volume_to_ref
    import numpy as np
    from torch.nn import Softmax

   # Set up threads and device
    device = torch.device(device)
    print('Using ' + args.device)
    if threads < 0:
        threads = os.cpu_count()
        print('Using all available threads ( %s )' % threads)
    else:
        print('Using %s thread(s)' % threads)
    torch.set_num_threads(threads)

    # Constants;  TODO:replace by FS paths
    model_file = os.path.join(os.environ.get('FREESURFER_HOME'), 'models', 'WMH-SynthSeg_v10_231110.pth')
    label_list_segmentation = [0, 14, 15, 16, 24, 77, 85, 2, 3, 4, 7, 8, 10, 11, 12, 13, 17, 18, 26, 28, 41, 42, 43, 46,
                               47, 49, 50, 51, 52, 53, 54, 58, 60]
    label_list_segmentation_torch = torch.tensor(label_list_segmentation, device=device)
    label_names = ['background', '3rd-ventricle', '4th-ventricle', 'brainstem', 'extracerebral_CSF', 'WMH', 'optic-chiasm',
                   'left-white-matter', 'left-cortex', 'left-lateral-ventricle', 'left-cerebellum-white-matter', 'left-cerebellum-cortex',
                   'left-thalamus', 'left-caudate', 'left-putamen', 'left-pallidum', 'left-hippocampus', 'left-amygdala', 'left-accumbens', 'left-ventral-DC',
                   'right-white-matter', 'right-cortex', 'right-lateral-ventricle', 'right-cerebellum-white-matter', 'right-cerebellum-cortex',
                   'right-thalamus', 'right-caudate', 'right-putamen', 'right-pallidum', 'right-hippocampus', 'right-amygdala', 'right-accumbens', 'right-ventral-DC']
    n_neutral_labels = 7
    n_labels = len(label_list_segmentation)
    in_channels = 1
    out_channels = n_labels + 4 + 1 + 1
    f_maps = 64
    layer_order = 'gcl'
    num_groups = 8
    num_levels = 5
    ref_res = np.array([1.0, 1.0, 1.0])
    voxelsize = np.prod(ref_res)

    # We enter PyTorch territory
    with torch.no_grad():

        # Model
        print('Preparing model and loading weights')

        model = UNet3D(in_channels, out_channels, final_sigmoid=False, f_maps=f_maps, layer_order=layer_order,
                       num_groups=num_groups, num_levels=num_levels, is_segmentation=False, is3d=True).to(device)
        checkpoint = torch.load(model_file)
        model.load_state_dict(checkpoint['model_state_dict'])

        n_ims = len(images_to_segment)
        if output_csv_path is not None:
            csv = open(output_csv_path, 'w')
            csv.write('Input-file,Intracranial-volume')
            for l in range(len(label_names)):
                lab = label_list_segmentation[l]
                if lab>0:
                    name = label_names[l]
                    csv.write(',' + name + '(' + str(lab) + ')')
            csv.write('\n')

        for nim in range(n_ims):
            input_file = images_to_segment[nim]
            output_file = segmentations_to_write[nim]
            print('Working on image ' + str(1+nim) + ' of ' + str(+n_ims) + ': ' + input_file)

            # try:

            print('     Loading input volume and normalizing to [0,1]')
            image, aff = MRIread(input_file)
            image_torch = torch.tensor(np.squeeze(image).astype(float), device=device)
            while len(image_torch.shape)>3:
                image_torch = image_torch.mean(image, dim=-1)
            image_torch, aff2 = align_volume_to_ref(image_torch, aff, aff_ref=np.eye(4), return_aff=True, n_dims=3)
            image_torch = image_torch / torch.max(image_torch)

            print('     Upscaling to target resolution')
            voxsize = np.sqrt(np.sum(aff2 ** 2, axis=0))[:-1]
            factors = voxsize / ref_res
            upscaled = myzoom_torch(image_torch, factors, device=device)
            aff_upscaled = aff2.copy()
            for j in range(3):
                aff_upscaled[:-1, j] = aff_upscaled[:-1, j] / factors[j]
            aff_upscaled[:-1, -1] = aff_upscaled[:-1, -1] - np.matmul(aff_upscaled[:-1, :-1], 0.5 * (factors - 1))
            upscaled_padded = torch.zeros(tuple((np.ceil(np.array(upscaled.shape) / 32.0) * 32).astype(int)), device=device)
            upscaled_padded[:upscaled.shape[0], :upscaled.shape[1], :upscaled.shape[2]] = upscaled

            print('     Pushing data through the CNN')
            pred1 = model(upscaled_padded[None, None, ...])[:, :, :upscaled.shape[0], :upscaled.shape[1], :upscaled.shape[2]]
            pred2 = torch.flip(model(torch.flip(upscaled_padded,[0])[None, None, ...]), [2])[:, :, :upscaled.shape[0], :upscaled.shape[1], :upscaled.shape[2]]
            softmax = Softmax(dim=0)
            nlat = int((n_labels - n_neutral_labels) / 2.0)
            vflip = np.concatenate([np.array(range(n_neutral_labels)),
                                    np.array(range(n_neutral_labels + nlat, n_labels)),
                                    np.array(range(n_neutral_labels, n_neutral_labels + nlat))])
            pred_seg_p = 0.5 * softmax(pred1[0, 0:n_labels, ...]) + 0.5 * softmax(pred2[0, vflip, ...])
            pred_seg = label_list_segmentation_torch[torch.argmax(pred_seg_p, 0)]
            pred_seg = np.squeeze(pred_seg.detach().cpu().numpy())

            # Write volumes from soft segmentation, if needed
            vols = voxelsize * torch.sum(pred_seg_p, dim=[1,2,3]).detach().cpu().numpy()
            if output_csv_path is not None:
                # Subject name and ICV
                csv.write(output_file + ',' + str(np.sum(vols[1:])))
                # volumes of structures
                for l in range(len(label_list_segmentation)):
                    if label_list_segmentation[l] > 0:
                        csv.write(',' + str(vols[l]) )
                csv.write('\n')

            print('     Writing segmentation to disk: ' + output_file)
            MRIwrite(pred_seg, aff_upscaled, output_file)
            print(' ')
            # except Exception as e:
            #     print("     An error occurred in this volume:", str(e))

        # We are done!
        if output_csv_path is not None:
            csv.close()
            print(' ')
            print('Written volumes to ' + output_csv_path)
            print(' ')
        print('All done!')
        print('If you use this method in a publication, please cite the following article:')
        print('Quantifying white matter hyperintensity and brain volumes in heterogeneous clinical and low-field portable MRI')
        print('Laso P, Cerri S, Sorby-Adams A, Guo J, Matteen F, Goebl P, Wu J, Li H, Young SI, Billot B, Puonti O, Rosen MS,')
        print('Kirsch J, Strisciuglio N, Wolterink JM, Eshaghi A, Barkhof F, Kimberly WT, and Iglesias JE')
        print('Under review.')


# execute script
if __name__ == '__main__':
    main()
