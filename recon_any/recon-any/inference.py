import sys
sys.path.append(__file__)
import torch
import os
from unet3d.model import UNet3D
import numpy as np
from torch.nn import Softmax
from argparse import ArgumentParser
from scipy.ndimage import binary_dilation
import nibabel as nib
from scipy.ndimage import gaussian_filter
from scipy.ndimage.morphology import distance_transform_edt
import interpol


# ================================================================================================
#                                         Main Entrypoint
# ================================================================================================

def main():

    # parse arguments
    parser = ArgumentParser()
    parser.add_argument("--input_image", type=str,
                        help="images to process. Can be the path to a single image or to a folder")
    parser.add_argument("--output_dir", type=str,
                        help="Directory where outputs will be written")
    parser.add_argument("--model_path", type=str,
                        help="path to model files")
    parser.add_argument("--case_type", type=str,
                        help="path to model files")
    parser.add_argument("--cpu", action="store_true", help="enforce running with CPU rather than GPU.")
    parser.add_argument("--threads", type=int, default=-1, dest="threads",
                        help="number of threads to be used by PyTorch when running on CPU (-1 for maximum).")


    args = vars(parser.parse_args())

    # TODO: change this at some point
    bspline_zooming = True

    # enforce CPU processing if necessary
    if args['cpu']:
        print('using CPU, hiding all CUDA_VISIBLE_DEVICES')
        os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
        device = 'cpu'
    else:
        device = 'cuda'

    # limit the number of threads to be used if running on CPU
    if args['threads'] < 0:
        args['threads'] = os.cpu_count()
        print('using all available threads ( %s )' % args['threads'])
    else:
        print('using %s thread(s)' % args['threads'])
    torch.set_num_threads(args['threads'])

    # Constants and models that depend on case type
    names_list = [None] * 100; names_list[0] = 'background'; names_list[14] = 'third-ventricle'; names_list[15] = 'fourth-ventricle';  names_list[16] = 'brainstem';
    names_list[24] = 'extracerebral-csf'; names_list[77] = 'wm-lesion'; names_list[85] = 'optic-chiasm'; names_list[2] = 'left-white-matter'; names_list[3] = 'left-cortex';
    names_list[4] = 'left-lateral-ventricle';  names_list[7] = 'left-cerebellum-white-matter';  names_list[8] = 'left-cerebellum-cortex';  names_list[10] = 'left-thalamus';
    names_list[11] = 'left-caudate';  names_list[12] = 'left-putamen';  names_list[13] = 'left-pallidum';  names_list[17] = 'left-hippocampus';  names_list[18] = 'left-amygdala';
    names_list[26] = 'left-accumbens';  names_list[28] = 'left-ventral-dc';  names_list[41] = 'right-white-matter';  names_list[42] = 'right-cortex';  names_list[43] = 'right-lateral-ventricle';
    names_list[46] = 'right-cerebellum-white-matter';  names_list[47] = 'right-cerebellum-cortex';  names_list[49] = 'right-thalamus';  names_list[50] = 'right-caudate';  names_list[51] = 'right-putamen';
    names_list[52] = 'right-pallidum';  names_list[53] = 'right-hippocampus';  names_list[54] = 'right-amygdala';  names_list[58] = 'right-accumbens';  names_list[60] = 'right-ventral-dc';
    label_list_segmentation = [0, 14, 15, 16, 24, 77, 85, 2, 3, 4, 7, 8, 10, 11, 12, 13, 17, 18, 26, 28, 41, 42, 43, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60]

    if args['case_type'] == 'left-c':
        label_list_segmentation = [0, 2, 3, 4, 10, 11, 12, 13, 17, 18, 26, 28, 77]
        n_neutral_labels = len(label_list_segmentation)
        model_file = os.path.join(args['model_path'], 'hemi.pth')
        fliplr = False
        n_labels = len(label_list_segmentation)
        out_channels = n_labels + 4
    elif args['case_type'] == 'left-ccb':
        label_list_segmentation = [0, 2, 3, 4, 7, 8, 10, 11, 12, 13, 16, 17, 18, 26, 28, 77]
        n_neutral_labels = len(label_list_segmentation)
        model_file = os.path.join(args['model_path'], 'hemi_with_cerebellum_and_brainstem.pth')
        fliplr = False
        n_labels = len(label_list_segmentation)
        out_channels = n_labels + 4
    elif args['case_type'] == 'right-c':
        label_list_segmentation = [0, 41, 42, 43, 49, 50, 51, 52, 53, 54, 58, 60, 77]
        n_neutral_labels = len(label_list_segmentation)
        model_file = os.path.join(args['model_path'], 'hemi.pth')
        fliplr = True
        n_labels = len(label_list_segmentation)
        out_channels = n_labels + 4
    elif args['case_type'] == 'right-ccb':
        label_list_segmentation = [0, 41, 42, 43, 46, 47, 49, 50, 51, 52, 16, 53, 54, 58, 60, 77]
        n_neutral_labels = len(label_list_segmentation)
        model_file = os.path.join(args['model_path'], 'hemi_with_cerebellum_and_brainstem.pth')
        fliplr = True
        n_labels = len(label_list_segmentation)
        out_channels = n_labels + 4
    elif args['case_type'] == 'both':
        label_list_segmentation = [0, 14, 15, 16, 24, 77, 85, 2, 3, 4, 7, 8, 10, 11, 12, 13, 17, 18, 26, 28, 41,
                                   42, 43, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60]
        n_neutral_labels = 7
        model_file = os.path.join(args['model_path'], 'full.pth')
        fliplr = False
        n_labels = len(label_list_segmentation)
        out_channels = n_labels + 6
    else:
        raise Exception('case_type must be left-c, left-ccb, right-c, right-ccb, or both')

    # some other constants
    max_surf_distance = 2.0
    nlat = vflip = dflip = None
    if args['case_type'] == 'both':
        nlat = int((n_labels - n_neutral_labels) / 2.0)
        vflip = np.concatenate([np.array(range(n_neutral_labels)),
                                np.array(range(n_neutral_labels + nlat, n_labels)),
                                np.array(range(n_neutral_labels, n_neutral_labels + nlat))])
        dflip = [n_labels+2, n_labels+3, n_labels, n_labels+1]


    # down to business
    with torch.no_grad():
        print('Reading, resampling, and padding input image')
        im, aff = MRIread(args['input_image'], im_only=False, dtype='float')
        if fliplr:
            aff[0, :] = -aff[0, :]
        im = torch.tensor(np.squeeze(im), dtype=torch.float32, device=device)
        if bspline_zooming: # TODO: make sure this is working
            factors = np.sqrt(np.sum(aff**2, axis=0))[:-1]
            im = interpol.resize(im, factor=list(factors), anchor='first', interpolation=3, bound='dct2', prefilter=True)
            for j in range(3):
                aff[:-1, j] = aff[:-1, j] / factors[j]
        else:
            im, aff = torch_resize(im, aff, 1.0, device)
        # import pdb; pdb.set_trace()
        im, aff = align_volume_to_ref(im, aff, aff_ref=np.eye(4), return_aff=True, n_dims=3)
        # im_orig = im.clone()
        while len(im.shape) > 3:  # in case it's rgb
            im = torch.mean(im, axis=-1)
        im = im - torch.min(im)
        im = im / torch.max(im)
        W = (np.ceil(np.array(im.shape) / 32.0) * 32).astype('int')
        idx = np.floor((W - im.shape) / 2).astype('int')
        S = torch.zeros(*W, dtype=torch.float32, device=device)
        S[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]] = im

        print('Preparing model and loading weights')
        model = UNet3D(1, out_channels=out_channels, final_sigmoid=False, f_maps=64, layer_order='gcl',
                               num_groups=8, num_levels=5, is_segmentation=False, is3d=True).to(device)
        if device=='cpu':
            checkpoint = torch.load(model_file, map_location = torch.device('cpu'))
        else:
            checkpoint = torch.load(model_file)
        model.load_state_dict(checkpoint['model_state_dict'])

        #############
        print('Pushing data through the CNN')
        if args['case_type'] == 'both':
            for p in range(2):
                if p==0:
                    print('  First pass: without flipping')
                else:
                    print('  Second pass: with flipping')
                    S = torch.flip(S, [0])
                pred = model(S[None, None, ...])
                if p==1:
                    pred = torch.flip(pred, [2])
                pred = pred[0, :, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
                softmax = Softmax(dim=0)
                if p==0:
                    pred_seg_p = 0.5 * softmax(pred[0:n_labels, ...])
                    pred_dist = 0.5 * pred[n_labels:n_labels+4, ...]
                    pred_im = 0.5 * pred[-2, ...]
                    pred_bf = 0.5 * pred[-1, ...]
                else:
                    pred_seg_p += 0.5 * softmax(pred[vflip, ...])
                    pred_dist += 0.5 * pred[dflip, ...]
                    pred_im += 0.5 * pred[-2, ...]
                    pred_bf += 0.5 * pred[-1, ...]
        else:

            pred = model(S[None, None, ...])
            pred = pred[0, :, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
            softmax = Softmax(dim=0)
            pred_seg_p = softmax(pred[0:n_labels, ...])
            pred_dist = torch.clamp(pred[n_labels:n_labels + 2, ...], min=-max_surf_distance, max=max_surf_distance)
            pred_im = pred[-2, ...]
            pred_bf = pred[-1, ...]

        # Postprocess
        pred_seg = torch.tensor(label_list_segmentation, device=device)[torch.argmax(pred_seg_p, 0)]
        pred_dist = torch.clamp(pred_dist, min=-max_surf_distance, max=max_surf_distance)
        # if len(im_orig.shape)==3:
        #     bf_corr = im_orig * torch.exp(-pred_bf)
        # else:
        #     bf_corr = im_orig * torch.exp(-pred_bf[..., None])


        # Detach from GPU
        # im_orig = im_orig.detach().cpu().numpy()
        pred_seg = pred_seg.detach().cpu().numpy()
        pred_dist = np.moveaxis(pred_dist.detach().cpu().numpy(), 0, -1)
        pred_im = pred_im.detach().cpu().numpy()
        # pred_bf = pred_bf.detach().cpu().numpy()
        # bf_corr = bf_corr.detach().cpu().numpy()

        # Make fake image with cortex
        a_dist_transform = 2
        b_dist_transform = 0.05
        W = pred_dist[:, :, :, 0]
        P = pred_dist[:, :, :, 1]
        fake =  10 + 35 * (1 - (np.tanh(a_dist_transform * (W+b_dist_transform)) + 1) / 2) + 65 * (1 - (np.tanh(a_dist_transform * (P+b_dist_transform)) + 1) / 2)
        if pred_dist.shape[3]>2:
            W = pred_dist[:, :, :, 2]
            P = pred_dist[:, :, :, 3]
            fake = fake + 10 + 35 * (1 - (np.tanh(a_dist_transform * (W+b_dist_transform)) + 1) / 2) + 65 * (1 - (np.tanh(a_dist_transform * (P+b_dist_transform)) + 1) / 2)

        # Blend with SynthSR for enhancement
        ssr = np.clip(pred_im * 110, 0, 120)
        mask_hp = ((pred_seg == 17) | (pred_seg == 53))
        mask_ct = ((pred_seg == 3) | (pred_seg == 42))
        mask = binary_dilation(mask_ct, iterations=2) & (~mask_hp)
        mindist = np.min(np.abs(pred_dist), axis=-1, keepdims=False)
        w = np.exp(-mindist[mask])
        norm = ssr.copy()
        norm[mask] = w * fake[mask] + (1-w) * ssr[mask]

        # Comput wm segmentation files
        WM = pred_seg.copy()
        WM[WM == 3] = 0
        WM[WM == 42] = 0
        WM[WM == 8] = 0
        WM[WM == 47] = 0
        WM[WM == 24] = 0
        WM[WM > 0] = 110
        WM[(pred_seg == 4) | (pred_seg == 43)] = 250

        FILLED = WM.copy()
        FILLED[pred_seg == 16] = 0
        FILLED[pred_seg == 7] = 0
        FILLED[pred_seg == 46] = 0
        if (args['case_type'] == 'right-ccb' or args['case_type'] == 'right-c'):
            FILLED[FILLED > 0] = 127
        elif (args['case_type'] == 'left-ccb' or args['case_type'] == 'left-c'):
            FILLED[FILLED > 0] = 255
        else: # both
            L = ((pred_seg > 0) & (pred_seg < 30) & (pred_seg < 75))
            R = ((pred_seg > 0) & (pred_seg > 30) & (pred_seg < 75))
            Dleft = distance_transform_edt(L == 0)
            Dright = distance_transform_edt(R == 0)
            FILLED[FILLED > 0] = 127
            FILLED[(FILLED > 0) & (Dleft < Dright)] = 255

        # Crucial to flip back if needed
        if fliplr:
            aff[0, :] = -aff[0, :]

        # OK in an ideal world this wouldn't be necessary but FS is just happier with this header

        refaff = np.array([[-1.0, 0.0, 0.0, 0.0],[0.0, 0.0, 1.0, 0.0],[0.0, -1.0, 0.0, 0.0],[0.0, 0.0, 0.0, 1.0]])
        pred_seg, _ = align_volume_to_ref(pred_seg, aff, aff_ref=refaff, return_aff=True, n_dims=3)
        norm, _ = align_volume_to_ref(norm, aff, aff_ref=refaff, return_aff=True, n_dims=3)
        ssr, _ = align_volume_to_ref(ssr, aff, aff_ref=refaff, return_aff=True, n_dims=3)
        WM, _ = align_volume_to_ref(WM, aff, aff_ref=refaff, return_aff=True, n_dims=3)
        FILLED, aff = align_volume_to_ref(FILLED, aff, aff_ref=refaff, return_aff=True, n_dims=3)

        # Registration
        print('Computing Talairach transform')
        estimate_and_write_talairach(pred_seg, aff, os.path.join(args['output_dir'], 'transforms', 'talairach.xfm'))

        # Write to disk and we're done!
        print('Writing to disk')
        if os.path.isdir(args['output_dir']) is False:
            os.mkdir(args['output_dir'])


        MRIwrite(pred_seg, aff, os.path.join(args['output_dir'], 'aseg.auto_noCCseg.mgz'))
        MRIwrite(norm, aff, os.path.join(args['output_dir'], 'norm.mgz'))
        MRIwrite(norm, aff, os.path.join(args['output_dir'], 'brain.mgz'))
        MRIwrite(norm, aff, os.path.join(args['output_dir'], 'brainmask.mgz'))
        MRIwrite(ssr, aff, os.path.join(args['output_dir'], 'SynthSR.mgz'))
        MRIwrite(WM, aff, os.path.join(args['output_dir'], 'wm.seg.mgz'))
        MRIwrite(FILLED, aff, os.path.join(args['output_dir'], 'filled.mgz'))

        # write teh volumes to disk
        vols = torch.sum(pred_seg_p, dim=[1, 2, 3]).detach().cpu().numpy() * np.abs(np.linalg.det(aff))
        f = open(os.path.join(args['output_dir'], '..','stats','SynthSeg.vols.csv'), 'w')
        f.write('structure-name, structure-label, volume_in_cubic_mm \n')
        for i in range(len(label_list_segmentation)-1): # skip background
            lab = label_list_segmentation[i+1]
            f.write(str(lab) + ', ' + names_list[lab] + ', ' +  str(vols[i+1]) + '\n')
        # f.write('total_intracranial_volume, ___, ' +  str(np.sum(vols[1:]))  +  '\n')
        f.close()

        print('All done')

# ================================================================================================
#                                         Auxiliary functions
# ================================================================================================

# Auxiliary function to fit affine matrices
def getM(ref, mov):

   zmat = np.zeros(ref.shape[::-1])
   zcol = np.zeros([ref.shape[1], 1])
   ocol = np.ones([ref.shape[1], 1])
   zero = np.zeros(zmat.shape)

   A = np.concatenate([
       np.concatenate([np.transpose(ref), zero, zero, ocol, zcol, zcol], axis=1),
       np.concatenate([zero, np.transpose(ref), zero, zcol, ocol, zcol], axis=1),
       np.concatenate([zero, zero, np.transpose(ref), zcol, zcol, ocol], axis=1)], axis=0)

   b = np.concatenate([np.transpose(mov[0, :]), np.transpose(mov[1, :]), np.transpose(mov[2, :])], axis=0)

   x = np.matmul(np.linalg.inv(np.matmul(np.transpose(A), A)), np.matmul(np.transpose(A), b))

   M = np.stack([
       [x[0], x[1], x[2], x[9]],
       [x[3], x[4], x[5], x[10]],
       [x[6], x[7], x[8], x[11]],
       [0, 0, 0, 1]])

   return M


###############

def estimate_and_write_talairach(imS, affS, filename):

    # Note: labels and centers of gravity precomputed from:
    # /autofs/space/panamint_005/users/iglesias/data/MNItemplates/mni_icbm152_nlin_sym_09c/mni_icbm152_t1_tal_nlin_sym_09c.nii.gz
    # /autofs/space/panamint_005/users/iglesias/data/MNItemplates/mni_icbm152_nlin_sym_09c/mni_icbm152_t1_tal_nlin_sym_09c.synthseg2.nii.gz
    labels = np.array(
        [2, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 24, 26, 28, 41, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58,
         60, 1001, 1002, 1003, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019,
         1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 2001, 2002,
         2003, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021,
         2022, 2023, 2024, 2025, 2026, 2027, 2028, 2029, 2030, 2031, 2032, 2033, 2034, 2035]).astype(int)
    atlCOG = np.array([[-27., -13., -31., -17., -24., -11., -14., -27., -20., 0., 0.,
                        0., -26., -23., 0., -8., -9., 27., 14., 33., 17., 25.,
                        11., 14., 27., 20., 26., 23., 9., 9., -54., -4., -37.,
                        -6., -23., -35., -44., -53., -6., -30., -23., -13., -5., -60.,
                        -23., -6., -49., -44., -49., -11., -48., -4., -44., -8., -4.,
                        -34., -10., -25., -54., -57., -9., -28., -45., -37., 55., 4.,
                        37., 6., 22., 35., 45., 52., 6., 31., 23., 13., 5.,
                        60., 23., 6., 49., 43., 50., 10., 48., 4., 43., 8.,
                        4., 33., 10., 25., 54., 57., 8., 30., 46., 37.],
                       [-18., -18., -13., -54., -63., -18., 12., 3., -2., -7., -46.,
                        -31., -20., -4., -21., 10., -16., -18., -21., -15., -54., -63.,
                        -18., 12., 3., -3., -20., -4., 10., -16., -43., 21., 11.,
                        -79., -4., -41., -67., -35., -46., -90., 29., -67., 42., -23.,
                        -32., -28., 16., 43., 32., -79., -19., -21., -7., -59., 38.,
                        49., 29., -63., -9., -35., 67., 13., -21., 2., -42., 20.,
                        13., -81., -6., -41., -66., -34., -46., -90., 29., -68., 41.,
                        -24., -32., -27., 16., 44., 32., -82., -19., -20., -7., -59.,
                        39., 50., 30., -62., -8., -34., 67., 13., -20., 2.],
                       [19., 15., -15., -35., -38., 6., 10., -1., -2., -4., -34.,
                        -34., -16., -20., 8., -9., -10., 19., 14., -14., -35., -38.,
                        6., 10., -1., -3., -16., -20., -9., -10., 8., 28., 49.,
                        20., -34., -21., 31., -24., 22., 0., -19., -5., -17., -13.,
                        -17., 58., 14., -14., 2., 7., 47., 40., 47., 38., 1.,
                        20., 47., 53., -4., 33., -11., -37., 9., -2., 6., 28.,
                        48., 21., -33., -21., 31., -25., 23., -1., -20., -5., -17.,
                        -13., -17., 59., 13., -14., 2., 7., 46., 40., 47., 39.,
                        1., 19., 47., 54., -5., 34., -11., -38., 8., -3.]]).astype(float)
    nlab = len(labels)

    refCOG = np.zeros([4, nlab])
    ok = np.ones(nlab)
    for l in range(nlab):
        aux = np.where(imS == labels[l])
        if len(aux[0]) > 50:
            refCOG[0, l] = np.median(aux[0])
            refCOG[1, l] = np.median(aux[1])
            refCOG[2, l] = np.median(aux[2])
            refCOG[3, l] = 1
        else:
            ok[l] = 0
    refCOG = np.matmul(affS, refCOG)[:-1, :]  # in RAS
    TAL = getM(refCOG[:, ok > 0], atlCOG[:, ok > 0])
    with open(filename, 'w') as f:
        f.write('MNI Transform File')
        f.write('\n')
        f.write('% avi2talxfm')
        f.write('\n')
        f.write('\n')
        f.write('Transform_Type = Linear;')
        f.write('\n')
        f.write('Linear_Transform = ')
        f.write('\n')
        f.write(str(TAL[0, 0]) + ' ' + str(TAL[0, 1]) + ' ' + str(TAL[0, 2]) + ' ' + str(TAL[0, 3]))
        f.write('\n')
        f.write(str(TAL[1, 0]) + ' ' + str(TAL[1, 1]) + ' ' + str(TAL[1, 2]) + ' ' + str(TAL[1, 3]))
        f.write('\n')
        f.write(str(TAL[2, 0]) + ' ' + str(TAL[2, 1]) + ' ' + str(TAL[2, 2]) + ' ' + str(TAL[2, 3]) + ';')
        f.write('\n')

####################




####################

def MRIwrite(volume, aff, filename, dtype=None):

    if dtype is not None:
        volume = volume.astype(dtype=dtype)

    if aff is None:
        aff = np.eye(4)
    header = nib.Nifti1Header()
    nifty = nib.Nifti1Image(volume, aff, header)

    nib.save(nifty, filename)

###############################

def MRIread(filename, dtype=None, im_only=False):

    assert filename.endswith(('.nii', '.nii.gz', '.mgz')), 'Unknown data file: %s' % filename

    x = nib.load(filename)
    volume = x.get_fdata()
    aff = x.affine

    if dtype is not None:
        volume = volume.astype(dtype=dtype)

    if im_only:
        return volume
    else:
        return volume, aff

############################


def myzoom_torch_anisotropic(X, aff, newsize, device):

    if len(X.shape)==3:
        X = X[..., None]

    factors = np.array(newsize) / np.array(X.shape[:-1])
    delta = (1.0 - factors) / (2.0 * factors)

    vx = torch.arange(delta[0], delta[0] + newsize[0] / factors[0], 1 / factors[0], dtype=torch.float, device=device)[:newsize[0]]
    vy = torch.arange(delta[1], delta[1] + newsize[1] / factors[1], 1 / factors[1], dtype=torch.float, device=device)[:newsize[1]]
    vz = torch.arange(delta[2], delta[2] + newsize[2] / factors[2], 1 / factors[2], dtype=torch.float, device=device)[:newsize[2]]

    vx[vx < 0] = 0
    vy[vy < 0] = 0
    vz[vz < 0] = 0
    vx[vx > (X.shape[0]-1)] = (X.shape[0]-1)
    vy[vy > (X.shape[1] - 1)] = (X.shape[1] - 1)
    vz[vz > (X.shape[2] - 1)] = (X.shape[2] - 1)

    fx = torch.floor(vx).int()
    cx = fx + 1
    cx[cx > (X.shape[0]-1)] = (X.shape[0]-1)
    wcx = vx - fx
    wfx = 1 - wcx

    fy = torch.floor(vy).int()
    cy = fy + 1
    cy[cy > (X.shape[1]-1)] = (X.shape[1]-1)
    wcy = vy - fy
    wfy = 1 - wcy

    fz = torch.floor(vz).int()
    cz = fz + 1
    cz[cz > (X.shape[2]-1)] = (X.shape[2]-1)
    wcz = vz - fz
    wfz = 1 - wcz

    Y = torch.zeros([newsize[0], newsize[1], newsize[2], X.shape[3]], dtype=torch.float, device=device)

    dtype = X.dtype
    for channel in range(X.shape[3]):
        Xc = X[:,:,:,channel]

        tmp1 = torch.zeros([newsize[0], Xc.shape[1], Xc.shape[2]], dtype=dtype, device=device)
        for i in range(newsize[0]):
            tmp1[i, :, :] = wfx[i] * Xc[fx[i], :, :] +  wcx[i] * Xc[cx[i], :, :]
        tmp2 = torch.zeros([newsize[0], newsize[1], Xc.shape[2]], dtype=dtype, device=device)
        for j in range(newsize[1]):
            tmp2[:, j, :] = wfy[j] * tmp1[:, fy[j], :] +  wcy[j] * tmp1[:, cy[j], :]
        for k in range(newsize[2]):
            Y[:, :, k, channel] = wfz[k] * tmp2[:, :, fz[k]] +  wcz[k] * tmp2[:, :, cz[k]]

    if Y.shape[3] == 1:
        Y = Y[:,:,:, 0]

    if aff is not None:
        aff_new = aff.copy()
        for c in range(3):
            aff_new[:-1, c] = aff_new[:-1, c] / factors[c]
        aff_new[:-1, -1] = aff_new[:-1, -1] - aff[:-1, :-1] @ (0.5 - 0.5 / factors)
        return Y, aff_new
    else:
        return Y

###############################

def torch_resize(I, aff, resolution, device, power_factor_at_half_width=5, dtype=torch.float32, slow=False):

    if torch.is_grad_enabled():
        with torch.no_grad():
            return torch_resize(I, aff, resolution, device, power_factor_at_half_width, dtype, slow)

    slow = slow or (device == 'cpu')
    voxsize = np.sqrt(np.sum(aff[:-1, :-1] ** 2, axis=0))
    newsize = np.round(I.shape[0:3] * (voxsize / resolution)).astype(int)
    factors = np.array(I.shape[0:3]) / np.array(newsize)
    k = np.log(power_factor_at_half_width) / np.pi
    sigmas = k * factors
    sigmas[sigmas<=k] = 0  # TODO: we could maybe remove this line, to make sure we always smooth a bit?

    if len(I.shape) not in (3, 4):
        raise Exception('torch_resize works with 3D or 3D+label volumes')
    no_channels = len(I.shape) == 3
    if no_channels:
        I = I[:, :, :, None]
    if torch.is_tensor(I):
        I = I.permute([3, 0, 1, 2])
    else:
        I = I.transpose([3, 0, 1, 2])

    It_lowres = None
    for c in range(len(I)):
        It = torch.as_tensor(I[c], device=device, dtype=dtype)[None, None]
        # Smoothen if needed
        for d in range(3):
            It = It.permute([0, 1, 3, 4, 2])
            if sigmas[d]>0:
                sl = np.ceil(sigmas[d] * 2.5).astype(int)
                v = np.arange(-sl, sl + 1)
                gauss = np.exp((-(v / sigmas[d]) ** 2 / 2))
                kernel = gauss / np.sum(gauss)
                kernel = torch.tensor(kernel,  device=device, dtype=dtype)
                if slow:
                    It = conv_slow_fallback(It, kernel)
                else:
                    kernel = kernel[None, None, None, None, :]
                    It = torch.conv3d(It, kernel, bias=None, stride=1, padding=[0, 0, int((kernel.shape[-1] - 1) / 2)])


        It = torch.squeeze(It)
        It, aff2 = myzoom_torch_anisotropic(It, aff, newsize, device)
        It = It.detach()
        if torch.is_tensor(I):
            It = It.to(I.device)
        else:
            It = It.cpu().numpy()
        if len(I) == 1:
            It_lowres = It[None]
        else:
            if It_lowres is None:
                if torch.is_tensor(It):
                    It_lowres = It.new_empty([len(I), *It.shape])
                else:
                    It_lowres = np.empty_like(It, shape=[len(I), *It.shape])
            It_lowres[c] = It

        torch.cuda.empty_cache()

    if not no_channels:
        if torch.is_tensor(I):
            It_lowres = It_lowres.permute([1, 2, 3, 0])
        else:
            It_lowres = It_lowres.transpose([1, 2, 3, 0])
    else:
        It_lowres = It_lowres[0]

    return It_lowres, aff2


@torch.jit.script
def conv_slow_fallback(x, kernel):
    """1D Conv along the last dimension with padding"""
    y = torch.zeros_like(x)
    x = torch.nn.functional.pad(x, [(len(kernel) - 1) // 2]*2)
    x = x.unfold(-1, size=len(kernel), step=1)
    x = x.movedim(-1, 0)
    for i in range(len(kernel)):
        y = y.addcmul_(x[i], kernel[i])
    return y



#######


def align_volume_to_ref(volume, aff, aff_ref=None, return_aff=False, n_dims=3):
    """This function aligns a volume to a reference orientation (axis and direction) specified by an affine matrix.
    :param volume: a numpy array
    :param aff: affine matrix of the floating volume
    :param aff_ref: (optional) affine matrix of the target orientation. Default is identity matrix.
    :param return_aff: (optional) whether to return the affine matrix of the aligned volume
    :param n_dims: number of dimensions (excluding channels) of the volume corresponding to the provided affine matrix.
    :return: aligned volume, with corresponding affine matrix if return_aff is True.
    """

    # work on copy
    aff_flo = aff.copy()

    # default value for aff_ref
    if aff_ref is None:
        aff_ref = np.eye(4)

    # extract ras axes
    ras_axes_ref = get_ras_axes(aff_ref, n_dims=n_dims)
    ras_axes_flo = get_ras_axes(aff_flo, n_dims=n_dims)

    # align axes
    aff_flo[:, ras_axes_ref] = aff_flo[:, ras_axes_flo]
    for i in range(n_dims):
        if ras_axes_flo[i] != ras_axes_ref[i]:
            if type(volume)==np.ndarray:
                volume = np.swapaxes(volume, ras_axes_flo[i], ras_axes_ref[i])
            else:
                volume = torch.swapaxes(volume, ras_axes_flo[i], ras_axes_ref[i])
            swapped_axis_idx = np.where(ras_axes_flo == ras_axes_ref[i])
            ras_axes_flo[swapped_axis_idx], ras_axes_flo[i] = ras_axes_flo[i], ras_axes_flo[swapped_axis_idx]

    # align directions
    dot_products = np.sum(aff_flo[:3, :3] * aff_ref[:3, :3], axis=0)
    for i in range(n_dims):
        if dot_products[i] < 0:
            if type(volume) == np.ndarray:
                volume = np.flip(volume, axis=i)
            else:
                volume = torch.flip(volume, [i])
            aff_flo[:, i] = - aff_flo[:, i]
            aff_flo[:3, 3] = aff_flo[:3, 3] - aff_flo[:3, i] * (volume.shape[i] - 1)

    if return_aff:
        return volume, aff_flo
    else:
        return volume

##############

def get_ras_axes(aff, n_dims=3):
    """This function finds the RAS axes corresponding to each dimension of a volume, based on its affine matrix.
    :param aff: affine matrix Can be a 2d numpy array of size n_dims*n_dims, n_dims+1*n_dims+1, or n_dims*n_dims+1.
    :param n_dims: number of dimensions (excluding channels) of the volume corresponding to the provided affine matrix.
    :return: two numpy 1d arrays of lengtn n_dims, one with the axes corresponding to RAS orientations,
    and one with their corresponding direction.
    """
    aff_inverted = np.linalg.inv(aff)
    img_ras_axes = np.argmax(np.absolute(aff_inverted[0:n_dims, 0:n_dims]), axis=0)
    return img_ras_axes


# execute script
if __name__ == '__main__':
    main()
