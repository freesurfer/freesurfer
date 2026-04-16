import sys
import os
sys.path.insert(0, os.path.join(os.environ.get('FREESURFER_HOME'),'python/packages/SuperSynth/'))
from argparse import ArgumentParser

# ================================================================================================
#                                         Main Entrypoint
# ================================================================================================

def main():

    # parse first
    parser = ArgumentParser(description='SuperSynth')
    parser.add_argument("--i", required=True, help="(required) Image to analyze, or CSV file with input,output,mode triplets.")
    parser.add_argument("--o", help="Output directory (ignored if input is CSV file)")
    parser.add_argument("--device", help="Device (cpu, cuda); default is cuda if available otherwise cpu")
    parser.add_argument("--sharpen_synths", action="store_true", help="(optional) Sharpen synthetic 1mm isotropic T1/T2/FLAIR with unsharp masking")
    parser.add_argument("--model_file", required=True, help="(required) Checkpoint .pth file")
    parser.add_argument("--test_time_flipping", action="store_true", help="(optional) Flipping for test-time augmentation.")
    parser.add_argument("--threads", type=int, default=-1, help="(optional) Number of CPU cores to be used (-1 = all available; this is the default)")
    parser.add_argument("--mode", help="Segmentation mode: invivo / cerebrum / left-hemi / right-hemi / exvivo (ignored if input is CSV file)")
    args = parser.parse_args()

    device = args.device
    threads = args.threads
    flipping = args.test_time_flipping
    model_file = args.model_file
    input_file = args.i
    outputdir = args.o
    mode = args.mode
    sharpen_synths = args.sharpen_synths

    # Don't bother importing stuff if parser fails hehe
    import torch
    import numpy as np
    import nibabel as nib
    from SuperSynth.utils import MRIread, MRIwrite, torch_resize, align_volume_to_ref, get_largest_connected_component, get_label_lists_etc, gaussian_blur_3d
    from SuperSynth.frugal_models import frugal_models
    from torch.nn.functional import softmax
    from scipy.ndimage.morphology import binary_dilation
    import csv
    import time
    import threading

    # Set up CPU/GPU and threads
    if device is None:
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    if device=='cpu':
        os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
    if threads < 0:
        threads = os.cpu_count()
    torch.set_num_threads(threads)


    # some constants
    label_list_segmentation_whole_freesurfer, label_list_segmentation_exvivo_freesurfer, label_list_segmentation_cerebrum_freesurfer, \
    label_list_segmentation_hemi_freesurfer_left, label_list_segmentation_hemi_freesurfer_right, label_list_segmentation_whole, \
    label_list_segmentation_hemis, label_list_segmentation_exvivo, n_neutral_labels_whole, n_neutral_labels_hemis, n_neutral_labels_exvivo, \
    n_neutral_labels_cerebrum, n_labels_whole, n_labels_hemis, n_labels_exvivo, n_labels_cerebrum, vflip_invivo, vflip_exvivo, vflip_cerebrum, \
    list_to_kill_photo_whole = get_label_lists_etc()

    # See if input is file or list of files, and prepare list of files accordingly
    input_file_list = []
    outputdir_list = []
    mode_list = []
    is_single_file = True
    try:
        aux = nib.load(input_file)
        print('It seems like your input is an image file')
    except:
        print('It seems like your input is not an image; assuming CSV file with list of triplets (and ignoring --o/--mode flags)')
        is_single_file = False

    if is_single_file:
        if outputdir is None:
            raise Exception('In single file mode, output file must be provided')
        if mode is None:
            raise Exception('In single file mode, the input variable --mode must be provided')
        input_file_list.append(input_file)
        outputdir_list.append(outputdir)
        mode_list.append(mode)
    else:
        with open(input_file, 'r') as file:
            csv_reader = csv.reader(file, delimiter=',')
            for row in csv_reader:
                if len(row)>0: # skip empty rows
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
        with torch.autocast(device_type='cuda', dtype=torch.float16):

            # Some more variables that we put in the GPU
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
            MNIqcseg_file = script_directory = os.path.dirname(os.path.abspath(__file__)) + '/../atlas/atlas.qc_seg.nii.gz'
            nets = frugal_models(model_file, MNIqcseg_file, device)

            def async_write(outputdir, aff, im, seg_discrete, T1, T2, FLAIR, FIELD, ribbon):
                MRIwrite(im, aff, outputdir + '/input_resampled.mgz', dtype=np.uint8)
                MRIwrite(seg_discrete, aff, outputdir + '/segmentation.mgz', dtype=np.uint16)
                MRIwrite(T1, aff, outputdir + '/SynthT1.mgz', dtype=np.uint8)
                MRIwrite(T2, aff, outputdir + '/SynthT2.mgz', dtype=np.uint8)
                MRIwrite(FLAIR, aff, outputdir + '/SynthFLAIR.mgz', dtype=np.uint8)
                MRIwrite(FIELD, aff, outputdir + '/mni_deformation.mgz', dtype=np.float32)
                MRIwrite(ribbon, aff, outputdir + '/ribbon.mgz', dtype=np.uint8)

            writer_thread = None

            # Start the stopwatch
            global_start_time = time.time()

            # Loop over images
            for im_idx in range(len(input_file_list)):

                  try:

                    start_time = time.time()

                    input_file = input_file_list[im_idx]
                    outputdir = outputdir_list[im_idx]
                    mode = mode_list[im_idx]

                    print('Working on image ' + str(im_idx+1) + ' of ' + str(len(input_file_list)) + ': ' + input_file)
                    print('   Mode is: ' + mode)

                    print('   Reading, resampling, and padding input image')
                    im, aff = MRIread(input_file, im_only=False, dtype='float')
                    im = torch.tensor(np.squeeze(im), dtype=torch.float32, device=device)
                    while len(im.shape) > 3:
                        im = im.mean(dim=-1)
                    im[im.isnan()] = 0
                    im, aff = torch_resize(im, aff, 1.0, device)
                    im, aff = align_volume_to_ref(im, aff, aff_ref=np.eye(4), return_aff=True, n_dims=3)
                    im /= im.max()
                    W = (np.ceil(np.array(im.shape) / 16.0) * 16).astype('int')
                    idx = np.floor((W - im.shape) / 2).astype('int')
                    S = torch.zeros(*W, dtype=torch.float32, device=device)
                    S[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]] = im

                    print('   Pushing data through the CNN')
                    pred = nets.ssynth_inference(S[None, None, ...])
                    reg = torch.permute(pred['reg'][0, :, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1],
                                        idx[2]:idx[2] + im.shape[2]], [1, 2, 3, 0])
                    T1 = pred['T1'][0, 0, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    T2 = pred['T2'][0, 0, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    FLAIR = pred['FLAIR'][0, 0, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    LP = pred['LP'][0, 0, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    LW = pred['LW'][0, 0, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    RP = pred['RP'][0, 0, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    RW = pred['RW'][0, 0, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    if mode == 'cerebrum':
                        seg = pred['seg'][0, mask_photo_or_cerebrum_whole, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    elif mode == 'left-hemi':
                        seg = pred['seg'][0, v_left, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    elif mode == 'right-hemi':
                        seg = pred['seg'][0, v_right, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    elif mode == 'exvivo':
                        seg = pred['seg'][0, mask_exvivo_whole, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    elif mode == 'invivo':
                        seg = pred['seg'][0, :, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                    else:
                        raise Exception('mode not supported: ' + mode)
                    seg = softmax(seg, dim=0)  # segmentations are activations, at this point

                    if flipping:
                        print('   Pushing flipped data through the CNN')
                        S = torch.flip(S, [0])
                        pred = nets.ssynth_inference(S[None, None, ...])

                        aux = torch.flip(torch.permute(pred['reg'][0, ...], [1, 2, 3, 0]), [0])
                        aux[..., 0] = -aux[..., 0]
                        reg = 0.5 * reg + 0.5 * aux[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2], :]
                        T1 = 0.5 * T1 + 0.5 * torch.flip(pred['T1'][0, 0, ...], [0])[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                        T2 = 0.5 * T2 + 0.5 * torch.flip(pred['T2'][0, 0, ...], [0])[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                        FLAIR = 0.5 * FLAIR + 0.5 * torch.flip(pred['FLAIR'][0, 0, ...], [0])[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1],idx[2]:idx[2] + im.shape[2]]
                        LP = 0.5 * LP + 0.5 * torch.flip(pred['RP'][0, 0, ...], [0])[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                        LW = 0.5 * LW + 0.5 * torch.flip(pred['RW'][0, 0, ...], [0])[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                        RP = 0.5 * RP + 0.5 * torch.flip(pred['LP'][0, 0, ...], [0])[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                        RW = 0.5 * RW + 0.5 * torch.flip(pred['LW'][0, 0, ...], [0])[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                        activations = torch.flip(pred['seg'][0, ...], [1])[:, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
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
                        else:  # 'invivo':
                            activations = activations[vflip_invivo, ...]
                        seg = 0.5 * seg + 0.5 * softmax(activations, dim=0)

                    print('   Postprocessing segmentation')
                    # Discretize segmentations
                    if mode == 'cerebrum':
                        seg_discrete = torch.tensor(label_list_segmentation_whole_freesurfer, device=device)[mask_photo_or_cerebrum_whole][torch.argmax(seg, 0)]
                    elif mode == 'left-hemi':
                        seg_discrete = torch.tensor(label_list_segmentation_hemi_freesurfer_left, device=device)[torch.argmax(seg, 0)]
                    elif mode == 'right-hemi':
                        seg_discrete = torch.tensor(label_list_segmentation_hemi_freesurfer_right, device=device)[torch.argmax(seg, 0)]
                    elif mode == 'exvivo':
                        seg_discrete = torch.tensor(label_list_segmentation_exvivo_freesurfer, device=device)[torch.argmax(seg, 0)]
                    elif mode == 'invivo':
                        seg_discrete = torch.tensor(label_list_segmentation_whole_freesurfer, device=device)[torch.argmax(seg, 0)]
                    else:
                        raise Exception('mode not supported: ' + mode)

                    print('   Postprocessing cortical ribbon')
                    a = 2
                    max_surf_distance= 3.0
                    LW = torch.clamp(LW, min=-max_surf_distance, max=max_surf_distance)
                    RW = torch.clamp(RW, min=-max_surf_distance, max=max_surf_distance)
                    LP = torch.clamp(LP, min=-max_surf_distance, max=max_surf_distance)
                    RP = torch.clamp(RP, min=-max_surf_distance, max=max_surf_distance)
                    ribbonL = 70 * (1 - (torch.tanh(a * (LW + 0.3)) + 1) / 2) + 40 * (1 - (torch.tanh(a * LP) + 1) / 2)
                    ribbonR = 70 * (1 - (torch.tanh(a * (RW + 0.3)) + 1) / 2) + 40 * (1 - (torch.tanh(a * RP) + 1) / 2)
                    if mode == 'left-hemi':
                        ribbon = ribbonL
                    elif mode == 'right-hemi':
                        ribbon = ribbonR
                    else:
                        ribbon = torch.maximum(ribbonL, ribbonR)

                    # get masks for fiting deformations and postprocessing segmentations
                    M = (seg_discrete > 0) & (seg_discrete != 24) & (seg_discrete < 900)  # useful for later
                    M = get_largest_connected_component(M.detach().cpu().numpy())
                    Mdilated = binary_dilation(M, iterations=2)
                    M = torch.tensor(M, device=device, dtype=torch.bool)
                    Mdilated = torch.tensor(Mdilated, device=device, dtype=torch.bool)
                    if mode != 'invivo':
                        T1[~Mdilated] = 0
                        T2[~Mdilated] = 0
                        FLAIR[~Mdilated] = 0
                        seg_discrete[~M] = 0
                    ribbon[~Mdilated] = 0

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
                            if (lab > 0) and ((lab != 24) or (mode=='invivo')) and (lab < 900) and (lab != 99):
                                if FSlabelNames[lab] is None:
                                    row1.append('(' + str(lab) + ')')
                                else:
                                    row1.append(FSlabelNames[lab] + ' (' + str(lab) + ')')
                                row2.append(str(vols[l]))
                        writer.writerow(row1)
                        writer.writerow(row2)

                    if sharpen_synths:
                        print('   Postprocessing Synth images')
                        amount_usm = 0.75
                        sigma_usm = 1.5
                        T1 += ((T1 - gaussian_blur_3d(T1, sigma_usm * np.ones(3), device)) * amount_usm)
                        T2 += ((T2 - gaussian_blur_3d(T2, sigma_usm * np.ones(3), device)) * amount_usm)
                        FLAIR += ((FLAIR - gaussian_blur_3d(T2, sigma_usm * np.ones(3), device)) * amount_usm)

                    print('   Fitting atlas transforms')

                    # First, affine
                    ri = torch.arange(reg.shape[0], dtype=torch.float32, device=device); mu_ri = ri.mean(); ri -= mu_ri; ri /= 100
                    rj = torch.arange(reg.shape[1], dtype=torch.float32, device=device); mu_rj = rj.mean(); rj -= mu_rj; rj /= 100
                    rk = torch.arange(reg.shape[2], dtype=torch.float32, device=device); mu_rk = rk.mean(); rk -= mu_rk; rk /= 100
                    mi, mj, mk = torch.meshgrid(ri, rj, rk, indexing='ij')
                    B = torch.stack([mi[M], mj[M], mk[M], torch.ones_like(mk[M])], dim=1)
                    P = torch.linalg.pinv(B)
                    fit_x = P @ reg[:, :, :, 0][M]; fit_y = P @ reg[:, :, :, 1][M]; fit_z = P @ reg[:, :, :, 2][M]
                    A1 = torch.eye(4, dtype=torch.float32, device=device); A1[0, 0] = 100; A1[1, 1] = 100; A1[2, 2] = 100
                    A2 = torch.eye(4, dtype=torch.float32, device=device); A2[0, -1] = mu_ri; A2[1, -1] = mu_rj; A2[2, -1] = mu_rk
                    A3 = torch.tensor(aff, dtype=torch.float32, device=device)
                    mni_affine = torch.stack([100 * fit_x, 100 * fit_y, 100 * fit_z, torch.tensor([0, 0, 0, 1], device=device)]) @ \
                                 (torch.linalg.inv(A1) @ torch.linalg.inv(A2) @ torch.linalg.inv(A3))
                    with open(outputdir + '/mni_affine.txt', 'w') as csvfile:
                        writer = csv.writer(csvfile)
                        for r in range(4):
                            row = []
                            for rr in range(4):
                                row.append(str(mni_affine[r, rr].detach().cpu().numpy()))
                            writer.writerow(row)

                    # Now, demons
                    print('   demons-like fit')
                    B = torch.stack([mi[Mdilated], mj[Mdilated], mk[Mdilated], torch.ones_like(mi[Mdilated])], dim=1)
                    aff_x = B @ fit_x; aff_y = B @ fit_y; aff_z = B @ fit_z
                    res_x = (reg[:, :, :, 0][Mdilated] - aff_x)
                    res_y = (reg[:, :, :, 1][Mdilated] - aff_y)
                    res_z = (reg[:, :, :, 2][Mdilated] - aff_z)
                    aux = torch.zeros(im.shape, dtype=res_x.dtype, device=device)
                    aux[:] = 0; aux[Mdilated] = res_x.clip(-.2, .2); res_x = 100 * gaussian_blur_3d(aux, [3, 3, 3], device)
                    aux[:] = 0; aux[Mdilated] = res_y.clip(-.2, .2); res_y = 100 * gaussian_blur_3d(aux, [3, 3, 3], device)
                    aux[:] = 0; aux[Mdilated] = res_z.clip(-.2, .2); res_z = 100 * gaussian_blur_3d(aux, [3, 3, 3], device)
                    aux[:] = 0; aux[Mdilated] = 1.0; aux = gaussian_blur_3d(aux, [3, 3, 3], device)
                    res_x /= aux; res_y /= aux; res_z /= aux  # nicely take care of edges of mask

                    print('   Benjamins QC network')
                    pred_qc = nets.qc_inference(seg_discrete, mode, aff, mni_affine)
                    if flipping:
                        pred_qc = 0.5 * pred_qc + 0.5 * nets.qc_inference(seg_discrete.flip([0]), mode, aff, mni_affine)
                    pred_qc_np = pred_qc.detach().cpu().numpy()
                    with open(outputdir + '/qc.csv', 'w') as csvfile:
                        writer = csv.writer(csvfile)
                        row1 = []
                        row2 = []
                        qc_names = ['general white matter', 'general grey matter', 'general csf', 'cerebellum',
                                    'brainstem', 'thalamus', 'putamen+pallidum', 'hippocampus+amygdala']
                        for l in range(len(pred_qc_np)):
                            row1.append(qc_names[l])
                            row2.append(str(pred_qc_np[l]))
                        writer.writerow(row1)
                        writer.writerow(row2)

                    print('   Writing to disk (asynchronously)')
                    im_np = (im * 255).clip(0, 255).detach().cpu().numpy()
                    seg_discrete_np = seg_discrete.detach().cpu().numpy()
                    T1_np = (100 * T1).clip(0, 255).detach().cpu().numpy()
                    T2_np = (100 * T2).clip(0, 255).detach().cpu().numpy()
                    FLAIR_np = (100 * FLAIR).clip(0, 255).detach().cpu().numpy()
                    FIELD_np = torch.stack([res_x, res_y, res_z], dim=-1).detach().cpu().numpy()
                    ribbon_np = ribbon.clip(0, 255).detach().cpu().numpy()

                    del im, seg_discrete, T1, T2, FLAIR, res_x, res_y, res_z, S, M, Mdilated, aux, mi, mj, mk, B, LP, RP, LW, RW, ribbon
                    if flipping:
                        del activations
                    torch.cuda.empty_cache()

                    if writer_thread is not None:
                        writer_thread.join()

                    # 4. start writing asynchronously
                    writer_thread = threading.Thread(
                        target=async_write,
                        args=(outputdir, aff, im_np, seg_discrete_np, T1_np, T2_np, FLAIR_np, FIELD_np, ribbon_np)
                    )
                    writer_thread.start()

                    print('   freeview ' + input_file + ' ' + outputdir + '/*.mgz & ')
                    print('   oocalc ' + outputdir + '/volumes.csv &')
                    print('   oocalc ' + outputdir + '/qc.csv &')
                    print('   cat ' + outputdir + '/mni_affine.txt')
                    print(f"   Processing this case took: {(time.time() - start_time):.4f} seconds")



                  except Exception as e:
                      print('*** Something went wrong with this file, skipping to the next ***')
                      print(f'Error trace: {e}')

            print('   All done, waiting for results of final image to get written do disk')
            if writer_thread:
                writer_thread.join()
            print(f"   Processing all cases took: {(time.time() - global_start_time):.4f} seconds")

#########################################################################################################3

# execute script
if __name__ == '__main__':
    main()

