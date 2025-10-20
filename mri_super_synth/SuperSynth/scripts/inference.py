import sys
import os
sys.path.insert(0, os.path.join(os.environ.get('FREESURFER_HOME'),'python/packages'))
from argparse import ArgumentParser
import torch
import numpy as np

# ================================================================================================
#                                         Main Entrypoint
# ================================================================================================

def main():

    # parse first
    parser = ArgumentParser(description='SuperSynth')
    parser.add_argument("--i", required=True, help="(required) Image to analyze, or CSV file with input,output,mode triplets.")
    parser.add_argument("--o", help="Output directory (ignored if input is CSV file)")
    parser.add_argument("--device", help="Device (cpu, cuda); default is cuda if available otherwise cpu")
    parser.add_argument("--force_tiling", action="store_true", help="(optional) Switch on tiled processing on the CPU")
    parser.add_argument("--model_file", required=True, help="(required) Checkpoint .pth file")
    parser.add_argument("--test_time_flipping", action="store_true", help="(optional) Flipping for test-time augmentation.")
    parser.add_argument("--threads", type=int, default=-1, help="(optional) Number of CPU cores to be used (-1 = all available; this is the default)")
    parser.add_argument("--mode", help="Segmentation mode: invivo / cerebrum / left-hemi / right-hemi / exvivo (ignored if input is CSV file)", required=True)
    args = parser.parse_args()

    device = args.device
    threads = args.threads
    flipping = args.test_time_flipping
    model_file = args.model_file
    input_file = args.i
    outputdir = args.o
    mode = args.mode
    force_tiling = args.force_tiling

    # Don't bother importing stuff if parser fails hehe
    from SuperSynth.ext.unet3d.model import EugeniosResidualEncoderUNet3D
    from SuperSynth.SuperSynth.utils import MRIread, MRIwrite, torch_resize, align_volume_to_ref
    from torch.nn import Softmax
    from SuperSynth.SuperSynth.generators import fast_3D_interp_torch
    from SuperSynth.SuperSynth.utils import get_largest_connected_component
    from scipy.ndimage.morphology import binary_dilation
    import csv
    import nibabel as nib

    # Set up CPU/GPU and threads
    if device is None:
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    if device=='cpu':
        os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
        tile_size = 160 if force_tiling else None
    else:
        tile_size = 160
    if threads < 0:
        threads = os.cpu_count()
    torch.set_num_threads(threads)


    # some constants
    mni_atlas_file = os.path.dirname(os.path.abspath(__file__)) + '/../atlas/MNI_atlas_sym.nii.gz'
    max_surf_distance = 3.0
    f_maps = 96
    label_list_segmentation_whole_freesurfer = [0, 14, 15, 16, 24, 77, 85, 99, 901, 902, 906, 907, 908, 909, 911, 912, 914, 915, 916,
                                                930, 2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 17, 18, 26, 819, 821, 843, 865, 869,
                                                41, 42, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58, 820, 822, 844, 866, 870]
    label_list_segmentation_exvivo_freesurfer = [0, 14, 15, 16, 77, 85, 99, 2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 17, 18, 26,
                                                 819, 821, 843, 865, 869, 41, 42, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58,
                                                 820, 822, 844, 866, 870]
    label_list_segmentation_cerebrum_freesurfer = [0,  77,  85,  99, 2, 3, 4, 5, 10, 11, 12, 13, 17, 18, 26, 819, 821, 843, 865, 869,
                                                   41,  42,  43,  44,  49,  50,  51,  52,  53,  54,  58, 820, 822, 844, 866, 870]
    label_list_segmentation_hemi_freesurfer_left = [0, 2, 3, 4, 5, 10, 11, 12, 13, 17, 18, 26, 77, 99, 819, 821, 843, 865, 869]
    label_list_segmentation_hemi_freesurfer_right = [0, 41, 42, 43, 44, 49, 50, 51, 52, 53, 54, 58, 77, 99, 820, 822, 844, 866, 870]
    label_list_segmentation_whole = [0, 11, 12, 13, 16, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46,
                                     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 15, 17, 47, 49, 51, 53, 55,
                                     18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 48, 50, 52, 54, 56]
    label_list_segmentation_hemis = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
    label_list_segmentation_exvivo = [0, 11, 12, 13, 31, 32, 33, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 15, 17, 34, 36, 38,
                                      40, 42, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 35, 37, 39, 41, 43]
    n_neutral_labels_whole = 20
    n_neutral_labels_hemis = len(label_list_segmentation_hemis)
    n_neutral_labels_exvivo = 7
    n_neutral_labels_cerebrum = 4
    n_labels_whole = len(label_list_segmentation_whole)
    n_labels_hemis = len(label_list_segmentation_hemis)
    n_labels_exvivo = len(label_list_segmentation_exvivo)
    n_labels_cerebrum = len(label_list_segmentation_cerebrum_freesurfer)
    nlat = int((n_labels_whole - n_neutral_labels_whole) / 2.0)
    vflip_invivo = np.concatenate([np.array(range(n_neutral_labels_whole)),
                                   np.array(range(n_neutral_labels_whole + nlat, n_labels_whole)),
                                   np.array(range(n_neutral_labels_whole, n_neutral_labels_whole + nlat))])

    nlat = int((len(label_list_segmentation_exvivo) - n_neutral_labels_exvivo) / 2.0)
    vflip_exvivo = np.concatenate([np.array(range(n_neutral_labels_exvivo)),
                                   np.array(range(n_neutral_labels_exvivo + nlat, len(label_list_segmentation_exvivo))),
                                   np.array(range(n_neutral_labels_exvivo, n_neutral_labels_exvivo + nlat))])
    nlat = int((len(label_list_segmentation_cerebrum_freesurfer) - n_neutral_labels_cerebrum) / 2.0)
    vflip_cerebrum = np.concatenate([np.array(range(n_neutral_labels_cerebrum)),
                                   np.array(range(n_neutral_labels_cerebrum + nlat, len(label_list_segmentation_cerebrum_freesurfer))),
                                   np.array(range(n_neutral_labels_cerebrum, n_neutral_labels_cerebrum + nlat))])
    final_layers = ['LP', 'LW', 'RP', 'RW', 'SR', 'BF', 'reg', 'seg', 'T1', 'T2', 'FLAIR', 'CT']
    final_layer_nf = [1, 1, 1, 1, 1, 1, 3, n_labels_whole, 1, 1, 1, 1]

    # See if input is file or list of files, and prepare list of files accordingly
    input_file_list = []
    outputdir_list = []
    mode_list = []
    try:
        aux = nib.load(input_file)
        print('It seems like your input is an image file')
        input_file_list.append(input_file)
        outputdir_list.append(outputdir)
        mode_list.append(mode)
    except:
        print('It seems like your input is not an image; assuming CSV file with list of triplets (and ignoring --o/--mode flags)')
        with open(input_file, 'r') as file:
            csv_reader = csv.reader(file, delimiter=',')
            for row in csv_reader:
                if row[0][0]!='#': # character to skip line
                    if os.path.exists(row[0]):
                        input_file_list.append(row[0])
                    else:
                        raise Exception('Input file does not exist: ' + row[0])
                    outputdir_list.append(row[1])
                    if (row[2]!='invivo') and (row[2]!='cerebrum') and (row[2]!='left-hemi') and (row[2]!='right-hemi') and (row[2]!='exvivo'):
                        print(row)
                        raise Exception('mode must be invivo/cerebrum/left-hemi/right-hemi/exvivo')
                    mode_list.append(row[2])

    # down to business
    with torch.no_grad():

        # Some more variables that we put in the GPU
        list_to_kill_photo_whole = [5, 6, 11, 12, 13, 16, 22, 23, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46]
        mask_photo_or_cerebrum_whole = torch.ones(len(label_list_segmentation_whole), dtype=torch.bool, device=device)
        for l in range(len(label_list_segmentation_whole)):
            if np.sum(np.array(list_to_kill_photo_whole) == label_list_segmentation_whole[l]) > 0:
                mask_photo_or_cerebrum_whole[l] = False
        v_left = []
        for lab in label_list_segmentation_hemi_freesurfer_left:
            v_left.append(np.where(np.array(label_list_segmentation_whole_freesurfer) == lab)[0][0])
        v_left = torch.tensor(v_left, device=device, dtype=torch.int32)
        v_right = []
        for lab in label_list_segmentation_hemi_freesurfer_right:
            v_right.append(np.where(np.array(label_list_segmentation_whole_freesurfer) == lab)[0][0])
        v_right = torch.tensor(v_right, device=device, dtype=torch.int32)
        mask_exvivo_whole = torch.ones(len(label_list_segmentation_whole), dtype=torch.bool, device=device)
        for l in range(len(label_list_segmentation_whole)):
            if np.sum(np.array(label_list_segmentation_exvivo_freesurfer) == label_list_segmentation_whole_freesurfer[l]) == 0:
                mask_exvivo_whole[l] = False

        # Load FreeSurfer labels
        d = {}
        with open(os.getenv('FREESURFER_HOME') + '/FreeSurferColorLUT.txt') as f:
            for l in f:
                if l.strip() and not l.lstrip().startswith('#'):
                    p = l.split()
                    if len(p) > 1 and p[0].isdigit(): d[int(p[0])] = p[1]
        FSlabelNames = [None] * (max(d) + 1) if d else []
        for k, v in d.items(): FSlabelNames[k] = v

        print('Preparing model and loading weights')
        if device == 'cpu':
            cp = torch.load(model_file, map_location=torch.device('cpu'))
        else:
            cp = torch.load(model_file)
        backbone = EugeniosResidualEncoderUNet3D(1, None, final_sigmoid=False, f_maps=f_maps, layer_order='cgl',
                        num_groups=8, num_levels=5, is_segmentation=False, is3d=True,
                        skip_final_convolution=True).to(device)
        backbone.load_state_dict(cp['backbone_state_dict'])
        final_convs = dict()
        for l in range(len(final_layers)):
            final_convs[final_layers[l]] = torch.nn.Conv3d(f_maps, final_layer_nf[l], 1, device=device)
            final_convs[final_layers[l]].load_state_dict(cp[final_layers[l] + '_state_dict'])

        # Loop over images
        for im_idx in range(len(input_file_list)):

              try:

                input_file = input_file_list[im_idx]
                outputdir = outputdir_list[im_idx]
                mode = mode_list[im_idx]

                print('Working on image ' + str(im_idx+1) + ' of ' + str(len(input_file_list)) + ': ' + input_file)
                print('   Mode is: ' + mode)

                print('   Reading, resampling, and padding input image')
                im, aff = MRIread(input_file, im_only=False, dtype='float')
                while len(im.shape)>3:
                    im = np.mean(im, axis=-1)
                im = torch.tensor(np.squeeze(im), dtype=torch.float32, device=device)
                im[im.isnan()] = 0
                im, aff = torch_resize(im, aff, 1.0, device)
                im, aff = align_volume_to_ref(im, aff, aff_ref=np.eye(4), return_aff=True, n_dims=3)
                im_orig = im.clone()
                while len(im.shape) > 3:  # in case it's rgb
                    im = torch.mean(im, axis=-1)
                im_maxi = torch.max(im)
                im = im  / im_maxi
                W = (np.ceil(np.array(im.shape) / 32.0) * 32).astype('int')
                if tile_size is not None:
                    W[W < tile_size] = tile_size
                idx = np.floor((W - im.shape) / 2).astype('int')
                S = torch.zeros(*W, dtype=torch.float32, device=device)
                S[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]] = im

                print('   Pushing data through the CNN')
                if (device=='cpu') and (force_tiling==False):
                    print('   Working on CPU; inference without tiling')
                    bb = backbone(S[None, None, ...])
                else:
                    print('   Working on ' + device + '; inference with tiling')
                    bb = process_tile(backbone, S[None, None, ...], tile_size=tile_size)

                if mode=='right-hemi':
                    LP = LW = None
                else:
                    LP = final_convs['LP'](bb)[0, 0, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    LW = final_convs['LW'](bb)[0, 0, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                if mode=='left-hemi':
                    RP = RW = None
                else:
                    RP = final_convs['RP'](bb)[0, 0, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    RW = final_convs['RW'](bb)[0, 0, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                reg = torch.permute(final_convs['reg'](bb)[0, :, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]], [1, 2, 3, 0])
                activations = final_convs['seg'](bb)[0, :, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                if mode=='cerebrum':
                    activations = activations[mask_photo_or_cerebrum_whole, ...]
                elif mode=='left-hemi':
                    activations = activations[v_left, ...]
                elif mode=='right-hemi':
                    activations = activations[v_right, ...]
                elif mode=='exvivo':
                    activations = activations[mask_exvivo_whole, ...]
                elif mode=='invivo':
                    pass
                else:
                    raise Exception('mode not supported: ' + mode)
                softmax = Softmax(dim=0)
                seg = softmax(activations)
                T1 = final_convs['T1'](bb)[0, 0, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                T2 = final_convs['T2'](bb)[0, 0, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                FLAIR = final_convs['FLAIR'](bb)[0, 0, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]

                if flipping:
                    S = torch.flip(S, [0])
                    if (device=='cpu') and (force_tiling==False):
                        bb = backbone(S[None, None, ...])
                    else:
                        bb = process_tile(backbone, S[None, None, ...], tile_size=tile_size)
                    if mode != 'right-hemi':
                        LP = 0.5 * LP + 0.5 * final_convs['RP'](bb)[0, 0].flip(dims=[0])[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                        LW = 0.5 * LW + 0.5 * final_convs['RW'](bb)[0, 0].flip(dims=[0])[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    if mode != 'left-hemi':
                        RP = 0.5 * RP + 0.5 * final_convs['LP'](bb)[0, 0].flip(dims=[0])[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                        RW = 0.5 * RW + 0.5 * final_convs['LW'](bb)[0, 0].flip(dims=[0])[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    aux = torch.flip(torch.permute(final_convs['reg'](bb)[0, ...], [1, 2, 3, 0]), [0])
                    aux[..., 0] = -aux[..., 0]
                    reg = 0.5 * reg + 0.5 * aux[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2], :]
                    activations = torch.flip(final_convs['seg'](bb)[0, ...], [1])[:, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    if mode == 'cerebrum':
                        activations = activations[mask_photo_or_cerebrum_whole, ...]
                        activations = activations[vflip_cerebrum, ...]
                    elif mode == 'left-hemi':
                        activations = activations[v_right, ...]
                    elif mode == 'right-hemi':
                        activations = activations[v_left, ...]
                    elif mode == 'exvivo':
                        activations = activations[mask_exvivo_whole, ...]
                        activations = activations[vflip_exvivo, ...]
                    else: # 'invivo':
                        activations = activations[vflip_invivo, ...]
                    seg = 0.5 * seg + 0.5 * softmax(activations)
                    T1 = 0.5 * T1 + 0.5 * torch.flip(final_convs['T1'](bb)[0, 0, ...], [0])[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    T2 = 0.5 * T2 + 0.5 * torch.flip(final_convs['T2'](bb)[0, 0, ...], [0])[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    FLAIR = 0.5 * FLAIR + 0.5 * torch.flip(final_convs['FLAIR'](bb)[0, 0, ...], [0])[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]

                # Basic postprocessing
                if LP is not None:
                    LP = torch.clamp(LP, min=-max_surf_distance, max=max_surf_distance)
                    LW = torch.clamp(LW, min=-max_surf_distance, max=max_surf_distance)
                if RP is not None:
                    RP = torch.clamp(RP, min=-max_surf_distance, max=max_surf_distance)
                    RW = torch.clamp(RW, min=-max_surf_distance, max=max_surf_distance)
                if mode=='cerebrum':
                    seg_discrete = torch.tensor(label_list_segmentation_whole_freesurfer, device=device)[mask_photo_or_cerebrum_whole][torch.argmax(seg, 0)]
                elif mode=='left-hemi':
                    seg_discrete = torch.tensor(label_list_segmentation_hemi_freesurfer_left, device=device)[torch.argmax(seg, 0)]
                elif mode=='right-hemi':
                    seg_discrete = torch.tensor(label_list_segmentation_hemi_freesurfer_right, device=device)[torch.argmax(seg, 0)]
                elif mode=='exvivo':
                    seg_discrete = torch.tensor(label_list_segmentation_exvivo_freesurfer, device=device)[torch.argmax(seg, 0)]
                elif mode=='invivo':
                    seg_discrete = torch.tensor(label_list_segmentation_whole_freesurfer, device=device)[torch.argmax(seg, 0)]
                else:
                    raise Exception('mode not supported: ' + mode)

                # get masks for fiting deformations and postprocessing segmentations
                M = (seg_discrete > 0) & (seg_discrete != 24) & (seg_discrete < 900) # useful for later
                M = get_largest_connected_component(M.detach().cpu().numpy())
                Mdilated = binary_dilation(M, iterations=2)
                M = torch.tensor(M, device=device, dtype=torch.bool)
                Mdilated = torch.tensor(Mdilated, device=device, dtype=torch.bool)
                if mode!='invivo':
                    T1[~Mdilated]=0
                    T2[~Mdilated]=0
                    FLAIR[~Mdilated]=0
                    seg_discrete[~M] = 0

                # postprocess soft segmentations and compute volumes
                seg[0][~Mdilated] = 1
                for l in range(seg.shape[0]):
                    seg[l][~Mdilated] = 0
                vols = seg.sum(dim=[1, 2, 3]).detach().cpu().numpy()
                if os.path.isdir(outputdir) is False:
                    os.mkdir(outputdir)
                with open(outputdir + '/volumes.csv', 'w') as csvfile:
                    writer = csv.writer(csvfile)
                    if mode == 'cerebrum':
                        llist = (torch.tensor(label_list_segmentation_whole_freesurfer, device=device)[mask_photo_or_cerebrum_whole]).detach().cpu().numpy()
                    elif mode == 'left-hemi':
                        llist = label_list_segmentation_hemi_freesurfer_left
                    elif mode == 'right-hemi':
                        llist = label_list_segmentation_hemi_freesurfer_right
                    elif mode == 'exvivo':
                        llist = label_list_segmentation_exvivo_freesurfer
                    elif mode == 'invivo':
                        llist = label_list_segmentation_whole_freesurfer
                    row1 = []
                    row2 = []
                    for l in range(len(llist)):
                        lab = llist[l]
                        if (lab > 0) and (lab != 24) and (lab < 900) and (lab != 99):
                            if FSlabelNames[lab] is None:
                                row1.append('(' + str(lab) + ')')
                            else:
                                row1.append(FSlabelNames[lab] + ' (' + str(lab) + ')')
                            row2.append(str(vols[l]))
                    writer.writerow(row1)
                    writer.writerow(row2)

                ########
                im_orig = im_orig.detach().cpu().numpy()
                seg_discrete = seg_discrete.detach().cpu().numpy()
                if LP is not None:
                    LP = LP.detach().cpu().numpy()
                    LW = LW.detach().cpu().numpy()
                if RP is not None:
                    RP = RP.detach().cpu().numpy()
                    RW = RW.detach().cpu().numpy()
                T1 = T1.detach().cpu().numpy()
                T2 = T2.detach().cpu().numpy()
                FLAIR = FLAIR.detach().cpu().numpy()
                if True:
                    print('   Deforming atlas')
                    MNI, aff2 = MRIread(mni_atlas_file)
                    A = np.linalg.inv(aff2)
                    MNI = torch.tensor(MNI, device=device, dtype=torch.float32)
                    A = torch.tensor(A, device=device, dtype=torch.float32)
                    xx = 100 * reg[:, :, :, 0][M]
                    yy = 100 * reg[:, :, :, 1][M]
                    zz = 100 * reg[:, :, :, 2][M]
                    ii = A[0, 0] * xx + A[0, 1] * yy + A[0, 2] * zz + A[0, 3]
                    jj = A[1, 0] * xx + A[1, 1] * yy + A[1, 2] * zz + A[1, 3]
                    kk = A[2, 0] * xx + A[2, 1] * yy + A[2, 2] * zz + A[2, 3]
                    vals = fast_3D_interp_torch(MNI, ii, jj, kk, 'linear', device)
                    DEF = torch.zeros_like(reg[..., 0])
                    DEF[M] = vals

                    print('   linear fit')
                    M[DEF==0] = False
                    ri = np.arange(reg.shape[0]).astype('float'); ri -= np.mean(ri); ri /= 100
                    rj = np.arange(reg.shape[1]).astype('float'); rj -= np.mean(rj); rj /= 100
                    rk = np.arange(reg.shape[2]).astype('float'); rk -= np.mean(rk); rk /= 100
                    mi, mj, mk = np.meshgrid(ri, rj, rk, sparse=False, indexing='ij')
                    mi = torch.tensor(mi, device=device, dtype=torch.float)[M]
                    mj = torch.tensor(mj, device=device, dtype=torch.float)[M]
                    mk = torch.tensor(mk, device=device, dtype=torch.float)[M]
                    B = torch.stack([mi, mj, mk, torch.ones_like(mk)], dim=1)
                    P = torch.linalg.pinv(B)
                    fit_x = P @ ii[vals>0]; fit_y = P @ jj[vals>0]; fit_z = P @ kk[vals>0]
                    iiaff = B @ fit_x; jjaff = B @ fit_y; kkaff = B @ fit_z
                    valsAff = fast_3D_interp_torch(MNI, iiaff, jjaff, kkaff, 'linear', device)
                    DEFaff = torch.zeros_like(reg[..., 0])
                    DEFaff[M] = valsAff

                    print('   Demons-like fit')
                    res_i = ii[vals > 0] - iiaff; res_j = jj[vals > 0] - jjaff;  res_k = kk[vals > 0] - kkaff
                    aux = torch.zeros_like(reg[..., 0]); aux[M] = res_i.clip(-20, 20); res_i = gaussian_blur_3d(aux, [3, 3, 3], device)
                    aux = torch.zeros_like(reg[..., 0]); aux[M] = res_j.clip(-20, 20); res_j = gaussian_blur_3d(aux, [3, 3, 3], device)
                    aux = torch.zeros_like(reg[..., 0]); aux[M] = res_k.clip(-20, 20); res_k = gaussian_blur_3d(aux, [3, 3, 3], device)
                    aux = torch.zeros_like(reg[..., 0]); aux[M] = 1.0; aux = gaussian_blur_3d(aux, [3, 3, 3], device)
                    res_i /= aux; res_j /= aux; res_k /= aux
                    valsDemons = fast_3D_interp_torch(MNI, iiaff+res_i[M], jjaff+res_j[M], kkaff+res_k[M], 'linear', device)
                    DEFdemons = torch.zeros_like(reg[..., 0])
                    DEFdemons[M] = valsDemons

                a = 2
                if LW is not None:
                    fakeL = 70 * (1 - (np.tanh(a * (LW + 0.3)) + 1) / 2) + 40 * (1 - (np.tanh(a * LP) + 1) / 2)
                else:
                    fakeL = np.zeros_like(RP)
                if RW is not None:
                    fakeR = 70 * (1 - (np.tanh(a * (RW + 0.3)) + 1) / 2) + 40 * (1 - (np.tanh(a * RP) + 1) / 2)
                else:
                    fakeR = np.zeros_like(LP)
                fake = fakeL + fakeR

                print('   Writing to disk')
                MRIwrite((im_orig/im_orig.max()*255).clip(0, 255), aff, outputdir + '/input_resampled.mgz', dtype=np.uint8)
                MRIwrite(seg_discrete, aff, outputdir + '/segmentation.mgz', dtype=np.uint16)
                MRIwrite((100 * T1).clip(0, 255), aff, outputdir + '/SynthT1.mgz', dtype=np.uint8)
                MRIwrite((100 * T2).clip(0, 255), aff, outputdir + '/SynthT2.mgz', dtype=np.uint8)
                MRIwrite((100 * FLAIR).clip(0, 255), aff, outputdir + '/SynthFLAIR.mgz', dtype=np.uint8)
                MRIwrite(fake.clip(0,255), aff, outputdir + '/fakeCortex.mgz', dtype=np.uint8)
                if True:
                    MRIwrite((100 * reg).detach().cpu().numpy(), aff, outputdir + '/mni_coordinates.mgz', dtype=np.float32)
                    MRIwrite(DEF.clip(0,255).detach().cpu().numpy(), aff, outputdir + '/mni_deformed_direct.mgz', dtype=np.uint8)
                    MRIwrite(DEFaff.clip(0,255).detach().cpu().numpy(), aff, outputdir + '/mni_deformed_affine.mgz', dtype=np.uint8)
                    MRIwrite(DEFdemons.clip(0,255).detach().cpu().numpy(), aff, outputdir + '/mni_deformed_demons.mgz', dtype=np.uint8)
                print('   freeview ' + input_file + ' ' + outputdir + '/*.mgz &')
                print('   oocalc ' + outputdir + '/volumes.csv &')
              except Exception as e:
                  print('*** Something went wrong with this file, skipping to the next ***')
                  print(f'Error trace: {e}')

        print('   All done')

#########################################################################################################3


def make_gaussian_kernel(sigma, device):
    sl = int(np.ceil(3 * sigma))
    ts = torch.linspace(-sl, sl, 2*sl+1, dtype=torch.float, device=device)
    gauss = torch.exp((-(ts / sigma)**2 / 2))
    kernel = gauss / gauss.sum()
    return kernel

def gaussian_blur_3d(input, stds, device):
    from torch.nn.functional import conv3d
    blurred = input[None, None, :, :, :]
    if stds[0]>0:
        kx = make_gaussian_kernel(stds[0], device=device)
        blurred = conv3d(blurred, kx[None, None, :, None, None], stride=1, padding=(len(kx) // 2, 0, 0))
    if stds[1]>0:
        ky = make_gaussian_kernel(stds[1], device=device)
        blurred = conv3d(blurred, ky[None, None, None, :, None], stride=1, padding=(0, len(ky) // 2, 0))
    if stds[2]>0:
        kz = make_gaussian_kernel(stds[2], device=device)
        blurred = conv3d(blurred, kz[None, None, None, None, :], stride=1, padding=(0, 0, len(kz) // 2))
    return torch.squeeze(blurred)


def process_tile(model, input_tensor, tile_size=160):
    B, C, H, W, D = input_tensor.shape
    weight_map = torch.zeros_like(input_tensor)
    # Blending mask
    x = y = z = torch.linspace(-1, 1, tile_size, device=input_tensor.device)
    xx, yy, zz = torch.meshgrid(x, y, z, indexing="ij")
    blending_mask = torch.cos(xx * torch.pi / 2) * torch.cos(yy * torch.pi / 2) * torch.cos(zz * torch.pi / 2)
    blending_mask = blending_mask[None, None, :, :].expand(1, C, -1, -1, -1)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                i1 = 0 if i == 0 else (H - tile_size)
                i2 = tile_size if i == 0 else H
                j1 = 0 if j == 0 else (W - tile_size)
                j2 = tile_size if j == 0 else W
                k1 = 0 if k == 0 else (D - tile_size)
                k2 = tile_size if k == 0 else D

                tile = input_tensor[:, :, i1:i2, j1:j2, k1:k2]
                with torch.no_grad():
                    tile_output = model(tile)
                if (i == 0) and (j == 0) and (k == 0):
                    output = torch.zeros([B, tile_output.shape[1], H, W, D], device=input_tensor.device)
                output[:, :, i1:i2, j1:j2, k1:k2] += tile_output * blending_mask
                weight_map[:, :, i1:i2, j1:j2, k1:k2] += blending_mask

    output /= torch.clamp(weight_map, min=1e-8)
    return output


# execute script
if __name__ == '__main__':
    main()
