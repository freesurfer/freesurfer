#Imports
import os
import sys
BASE_PATH = os.path.join(os.environ.get('FREESURFER_HOME'),'python/packages/ERC_bayesian_segmentation/')
sys.path.insert(0, BASE_PATH)
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
import ext.my_functions as my
from datetime import datetime
import ERC_bayesian_segmentation.relabeling as relab
import glob
import scipy.sparse as sp
import ext.bias_field_correction_torch as bf
import csv
import argparse
import math
from torch.nn.functional import grid_sample
from ext.fireants_trimmed import Image, BatchedImages
from ext.fireants_trimmed import GreedyRegistration
import SimpleITK as sitk
from time import time
from ext.fireants_trimmed import HybridDiceLabelDiffloss
from scipy.ndimage import binary_erosion, binary_fill_holes
from ERC_bayesian_segmentation import SuperSynth_inference

########################################################

parser = argparse.ArgumentParser(description='Bayesian segmentation.')
parser.add_argument("--i", help="Image to segment.")
parser.add_argument("--model_file", help="Multi-task deep learning model for preprocessing")
parser.add_argument("--atlas_dir", help="Atlas directory")
parser.add_argument("--mode", help="Type of input (invivo, exvivo, cerebrum, hemi).")
parser.add_argument("--side", help="Hemisphere to segment (left or right).")
parser.add_argument("--bf_mode", help="bias field basis function: dct, polynomial, or hybrid", default="dct")
parser.add_argument("--o", help="Output directory.")
parser.add_argument("--write_rgb", action="store_true", help="Write soft segmentation to dis as RGB file.")
parser.add_argument("--write_bias_corrected", action="store_true", help="Write bias field corrected image to disk")
parser.add_argument("--device", help="Device (cpu, cuda)")
parser.add_argument("--device_registration", help="Use this option if you want to use a different device for the registration (useful for high-res ex vivo)")
parser.add_argument("--threads", type=int, default=-1, help="(optional) Number of CPU cores to be used. Default is -1 (use all available cores")
parser.add_argument("--skip", type=int, default=1, help="(optional) Skipping factor to easy memory requirements of priors when estimating Gaussian parameters. Default is 1.")
parser.add_argument("--resolution", type=float, default=0.4, help="(optional) Resolution of output segmentation")
parser.add_argument("--force_tiling", action="store_true", help="Forces tiling on CPU so it gives the same result as GPU")
parser.add_argument("--skip_bf", action="store_true", help="Skip bias field correction")
parser.add_argument("--smooth_grad_sigma", type=float, default=1.00, help="(optional) Parameter of Greedy FireANTs registration")
parser.add_argument("--smooth_warp_sigma", type=float, default=0.25, help="(optional) Parameter of Greedy FireANTs registration")
parser.add_argument("--optimizer_lr", type=float, default=0.5, help="(optional) Parameter of Greedy FireANTs registration")
parser.add_argument("--cc_kernel_size", type=int, default=7, help="(optional) Parameter of Greedy FireANTs registration")
parser.add_argument("--rel_weight_labeldiff", type=float, default=2.5, help="(optional) Relative weight of labels in Greedy FireANTs registration")
parser.add_argument("--save_atlas_nonlinear_reg", action="store_true", help="Save nonlinear atlas registration")
args = parser.parse_args()

########################################################

if args.i is None:
    raise Exception('Input image is required')
if args.model_file is None:
    raise Exception('Model file is required')
if args.atlas_dir is None:
    raise Exception('Atlas directory must be provided')
if args.o is None:
    raise Exception('Output directory must be provided')
if args.side is None:
    raise Exception('side must be provided (left or right)')
else:
    if (args.side!='left') and (args.side!='right'):
        raise Exception('Side must be left or right, but you specified: ' + args.side)
if args.mode is None:
    raise Exception('mode must be provided (invivo/cerebrum/hemi/exvivo)')
else:
    if (args.mode!='invivo') and (args.mode!='cerebrum') and (args.mode!='hemi') and (args.mode!='exvivo'):
        raise Exception('Mode must be invivo/cerebrum/hemi/exvivo, but you specified: ' + args.mode)
if (args.resolution<=0):
    raise Exception('Resolution must be non-negative; Exitting...')
if (args.skip<1):
    raise Exception('Skip cannot be less than 1; exitting...')


########################################################

if torch.cuda.is_available():
    print('GPU/CUDA seems to be available')
    if args.device is None:
        args.device = 'cpu'
    if args.device_registration is None:
        args.device_registration = args.device
    print('You selected use of the following devices:')
    print('  Registration: ' + args.device_registration)
    device_registration = torch.device(args.device_registration)
    print('  Rest of code: ' + args.device)
    device = torch.device(args.device)
    os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "garbage_collection_threshold:0.5,max_split_size_mb:32"

else:
    print('GPU/CUDA not available; using the CPU')
    device = torch.device('cpu')
    device_registration = torch.device('cpu')

########################################################

# limit the number of CPU threads to be used
if args.threads<0:
    args.threads = os.cpu_count()
    print('using all available CPU threads ( %s )' % args.threads)
else:
    print('using %s CPU thread(s)' % args.threads)
torch.set_num_threads(args.threads)

########################################################

# Reproducibility: https://discuss.pytorch.org/t/reproducibility-with-all-the-bells-and-whistles/81097
seed = 0; torch.manual_seed(seed); torch.cuda.manual_seed_all(seed); torch.cuda.manual_seed(seed)
np.random.seed(seed); np.random.seed(seed)
torch.backends.cudnn.deterministic = True; torch.backends.cudnn.benchmark = False

########################################################

# Disable gradients
torch.set_grad_enabled(False)

################

# Input data
input_volume = args.i
model_file = args.model_file
atlas_dir = args.atlas_dir
LUT_file = os.path.join(BASE_PATH, 'data_simplified', 'AllenAtlasLUT')
output_dir = args.o
skip_bf = args.skip_bf
bf_mode = args.bf_mode
side = args.side
mode = args.mode
resolution = args.resolution
skip = args.skip
force_tiling =  args.force_tiling

########################################################
# Detect problems with output directory right off the bat
if os.path.exists(output_dir):
    if len(os.listdir(output_dir))==0:
        print('Warning: output directory exists')
    else:
        pass
        #raise Exception('Ouput directory exists and is not empty; exitting...')
else:
    os.mkdir(output_dir)


########################################################

# Constants
dtype = torch.float32
SET_BG_TO_CSF = True # True = median of ventricles -> it seems much better than 0!
RESOLUTION_ATLAS = 0.2
TOL = 1e-9

############

if dtype == torch.float64:
    numpy_dtype = np.float64
elif dtype == torch.float32:
    numpy_dtype = np.float32
elif dtype == torch.float16:
    numpy_dtype = np.float16
else:
    raise Exception('type not supported')

########################################################

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

########################################################
print('Reading input image')
Iim, aff = my.MRIread(input_volume)
Iim = np.squeeze(Iim)
while len(Iim.shape) > 3:
    Iim = np.mean(Iim, axis=-1)

########################################################
print('Analyzing image with neural network...')
print('  Resampling, reorienting, and padding')
im = torch.tensor(Iim, dtype=dtype, device=device).squeeze()
im, imaff = my.torch_resize(im, aff, 1.0, device)
im, imaff = my.align_volume_to_ref(im, imaff, aff_ref=np.eye(4), return_aff=True, n_dims=3)
im /= im.max()
mode_supersynth = (side + '-hemi') if mode=='hemi' else mode
seg, reg = SuperSynth_inference.run_inference(im, True, model_file, mode_supersynth, output_dir + '/supersynth.vols.csv', device, force_tiling=force_tiling)
my.MRIwrite(seg.detach().cpu().numpy(), imaff, output_dir + '/supersynth.nii.gz')

########################################################

# Kill contralateral labels if needed
if mode=='hemi':
    print('Mode is single hemi, no need to kill contralateral labels')
else:
    print('Killing labels in contralateral hemisphere')
    if side=='left':
        labels_to_kill =  [0, 14, 15, 41, 42, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58, 820, 822, 844, 866, 870]
    else:
        labels_to_kill =  [0, 14, 15,  2,  3,  4,  5,  7,  8, 10, 11, 12, 13, 17, 18, 26, 819, 821, 843, 865, 869]
    for l in labels_to_kill:
        seg[seg == l] = 0


########################################################

print('   Reslice segmentation and left-right map to the space of the input image')
II, JJ, KK = np.meshgrid(np.arange(Iim.shape[0]), np.arange(Iim.shape[1]), np.arange(Iim.shape[2]), indexing='ij')
II = torch.tensor(II, device=device, dtype=dtype)
JJ = torch.tensor(JJ, device=device, dtype=dtype)
KK = torch.tensor(KK, device=device, dtype=dtype)
affine = np.linalg.inv(imaff) @ aff
II2 = affine[0, 0] * II + affine[0, 1] * JJ + affine[0, 2] * KK + affine[0, 3]
JJ2 = affine[1, 0] * II + affine[1, 1] * JJ + affine[1, 2] * KK + affine[1, 3]
KK2 = affine[2, 0] * II + affine[2, 1] * JJ + affine[2, 2] * KK + affine[2, 3]
Sim = my.fast_3D_interp_torch(seg, II2, JJ2, KK2, 'nearest').detach().cpu().numpy()
LRmap = my.fast_3D_interp_torch(reg[...,0], II2, JJ2, KK2, 'linear').detach().cpu().numpy()
APmap = my.fast_3D_interp_torch(reg[...,1], II2, JJ2, KK2, 'linear').detach().cpu().numpy()
ISmap = my.fast_3D_interp_torch(reg[...,2], II2, JJ2, KK2, 'linear').detach().cpu().numpy()

########################################################

if skip_bf==False:
    print('Correcting bias field')
    print('   Trying model with polynomial basis functions')
    try:
        Iim, _, bfmask = bf.correct_bias(Iim, Sim, maxit=100, penalty=0.1, order=6, device=device, dtype=dtype, basis=bf_mode, dontmask=True, cerebrum_only=((mode=='hemi') or (mode=='cerebrum')))
    except:
        if device.type=='cpu':
            raise Exception('Bias correction failed (out of memory?)')
        else:
            print('Bias correction on GPU failed; trying with CPU')
            Iim, _, bfmask = bf.correct_bias(Iim, Sim, maxit=100, penalty=0.1, order=4, device='cpu', dtype=dtype, basis=bf_mode, dontmask=True, cerebrum_only=((mode=='hemi') or (mode=='cerebrum')))
    if args.write_bias_corrected:
        bfmask = binary_fill_holes(bfmask)
        aux = Iim.copy()
        aux[~bfmask] = 0
        aux = (aux / np.max(Iim[bfmask]) * 255).astype(np.uint8)
        my.MRIwrite(aux, aff, output_dir + '/bias.corrected.nii.gz')
        del aux, bfmask

print('Normalizing intensities')
Iim = Iim * 110 / np.median(Iim[(Sim==2) | (Sim==41)])

# We should do tensors at this point...
Sim = torch.tensor(Sim, dtype=torch.int, device=device)
Iim = torch.tensor(Iim, dtype=dtype, device=device)
LRmap = torch.tensor(LRmap, dtype=dtype, device=device)
APmap = torch.tensor(APmap, dtype=dtype, device=device)
ISmap = torch.tensor(ISmap, dtype=dtype, device=device)

########################################################
# Fit affine registration to atlas
print('Linear fit of predicted MNI coordinates')
Mfit = ((Sim>0) & ((Sim<900) | (Sim>1000))  & (Sim!=4)  & (Sim!=5)  & (Sim!=43)  & (Sim!=44)  & (Sim!=14)  & (Sim!=15)  & (Sim!=24))
Mfit = torch.tensor(binary_erosion(Mfit.detach().cpu().numpy(), iterations=2), device=device, dtype=torch.bool)
# avoid numerical issues, we don't need a trillion voxels to fit this
prop = 100000 / Mfit.sum()
if (prop<1):
    aux = (torch.rand(Mfit.shape, device=device) < prop)
    Mfit = (Mfit & aux)
    del aux
ri = np.arange(Sim.shape[0]).astype('float'); mu_ri = np.mean(ri); ri -= mu_ri ; ri /= 100
rj = np.arange(Sim.shape[1]).astype('float'); mu_rj = np.mean(rj); rj -= mu_rj; rj /= 100
rk = np.arange(Sim.shape[2]).astype('float'); mu_rk = np.mean(rk); rk -= mu_rk; rk /= 100
mi, mj, mk = np.meshgrid(ri, rj, rk, sparse=False, indexing='ij')
mi = torch.tensor(mi, device=device, dtype=dtype)[Mfit]
mj = torch.tensor(mj, device=device, dtype=dtype)[Mfit]
mk = torch.tensor(mk, device=device, dtype=dtype)[Mfit]
B = torch.stack([mi, mj, mk, torch.ones_like(mk)], dim=1)
P = torch.linalg.pinv(B)
fit_lr = P @ (100*LRmap[Mfit]); fit_ap = P @ (100*APmap[Mfit]); fit_is = P @ (100*ISmap[Mfit])
aux = torch.stack([fit_lr, fit_ap, fit_is, torch.tensor([0,0,0,1], device=device, dtype=dtype)])
mat1 = np.matrix('1 0 0 ' + str(-mu_ri) + '; 0 1 0 ' + str(-mu_rj) + '; 0 0 1 ' + str(-mu_rk) + '; 0 0 0 1')
mat2 = np.diag([0.01, 0.01, 0.01, 1])
mat3 = torch.stack([fit_lr, fit_ap, fit_is, torch.tensor([0,0,0,1], device=device, dtype=dtype)]).detach().cpu().numpy()
M_input_vox_to_mni_ras = mat3 @ mat2 @ mat1
# MNI, affmni2 = my.MRIread('/homes/2/iglesias/gca.mgz'); my.MRIwrite(MNI, aff @ np.linalg.inv(M_input_vox_to_mni_ras)  @ affmni2, '/tmp/test.mgz')

########################################################
# Kill bottom of medulla and subdivide brainstem if needed
if (mode=='invivo') or (mode=='exvivo'):
    print('Dealing with bilateral labels: subdividing brainstem, optic chiasm, lesions');
    print('  (we also crop the bottom of the brainstem a bit)')
    LEFT = (LRmap < 0)
    Sim[(Sim == 16) & (ISmap < (-0.60))] = 0
    if (side == 'left'):
        Sim[(Sim == 16) & LEFT] = 161
        Sim[(Sim == 16)] = 0
        Sim[(Sim == 77) & (~LEFT)] = 0
        Sim[(Sim == 85) & (~LEFT)] = 0
    else:
        Sim[(Sim == 16) & (~LEFT)] = 162
        Sim[(Sim == 16)] = 0
        Sim[(Sim == 77) & LEFT] = 0
        Sim[(Sim == 85) & LEFT] = 0
else:
    print('No need to deal with bilateral labels')

# Release memory
del LRmap, APmap, ISmap, reg, seg, II, JJ, KK, II2, JJ2, KK2, ri, rj, rk, mi, mj, mk, B, P

#######################################

# Prepare data for hemisphere at hand
print('  Creating mask for tissue to segment (leave out ventricles)')
M = (Sim>0)
M[Sim==4] = 0
M[Sim==5] = 0
M[Sim==43] = 0
M[Sim==44] = 0
M[Sim==14] = 0
M[Sim==15] = 0
M = torch.tensor(my.getLargestCC(M.detach().cpu().numpy()), device=M.device, dtype=torch.bool)

# I now do this with the resampled mask later on, to avoid blurring mask edges
Mim, cropping = my.cropLabelVolTorch(M, margin=5)
del M
Sim = my.applyCropping(Sim, cropping)
Iim = my.applyCropping(Iim, cropping)
Iim_before_masking = Iim.clone()
Iim[Mim==0] = 0
aff[:3, -1] = aff[:3, -1] + aff[:-1, :-1] @ np.array([cropping[0].detach().cpu().numpy(), cropping[1].detach().cpu().numpy(), cropping[2].detach().cpu().numpy()])


########################################################

# Read atlas
print('Reading in atlas')

# ASEG labels are hard coded for each hemi, no biggie
if (side=='left'):
    aseg_label_list = np.array(
        [2, 3, 7, 8, 10, 11, 12, 13, 17, 18, 26, 28, 77, 85, 161, 819, 821, 843, 865, 869]).astype(int)
else:
    aseg_label_list = np.array(
        [41, 42, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60, 77, 85, 162, 820, 822, 844, 866, 870]).astype(int)

# Get the label groupings and atlas labels from the config files
tissue_index, grouping_labels, label_list, number_of_gmm_components = relab.get_tissue_settings(
            os.path.join(BASE_PATH, 'data_simplified', 'atlas_names_and_labels.yaml'),
            os.path.join(BASE_PATH, 'data_simplified', 'combined_atlas_labels_fireants.yaml'),
            os.path.join(BASE_PATH, 'data_simplified', 'combined_aseg_labels_new_targets.yaml'),
            os.path.join(BASE_PATH, 'data_simplified', 'gmm_components_fireants.yaml'),
            aseg_label_list
)

tidx = tissue_index[np.where(label_list == 0)[0][0]]
if tidx>0:
    raise Exception('First tissue class must be the background')
n_tissues = np.max(tissue_index) + 1
n_labels = len(label_list)
atlas_names = sorted(glob.glob(atlas_dir + '/label_*.npz'))
atlas_size = np.load(atlas_dir + '/size.npy')

class LabelDataset(Dataset):

    def __init__(self, fnames):
        self.fnames = fnames

    def __len__(self):
        return len(self.fnames)

    def __getitem__(self, item):
        print(item, self.fnames[item])
        prior = sp.load_npz(self.fnames[item])
        prior_indices = torch.as_tensor(prior.row)
        prior_values = torch.as_tensor(prior.data)
        return prior_indices, prior_values

# TODO: without this line, I get weird runtime errors...
#Haiden's Edit:
prefetch = 1
workers = 0
label_loader = DataLoader(LabelDataset(atlas_names), num_workers=workers)
A = np.zeros([*atlas_size, n_tissues], dtype=numpy_dtype)
# We keep track of these probability masses we use to correct differences in labeling between A and SynthSeg
MLhippo = np.zeros(atlas_size, dtype=numpy_dtype)
CLAUSTRUM = np.zeros(atlas_size, dtype=numpy_dtype)
RETICULAR = np.zeros(atlas_size, dtype=numpy_dtype)
LGN =  np.zeros(atlas_size, dtype=numpy_dtype)
AMYGDALA =  np.zeros(atlas_size, dtype=numpy_dtype)
FORNIX =  np.zeros(atlas_size, dtype=numpy_dtype)
CHIASM =  np.zeros(atlas_size, dtype=numpy_dtype)
MAMBODY =  np.zeros(atlas_size, dtype=numpy_dtype)
SEPTAL =  np.zeros(atlas_size, dtype=numpy_dtype)
for n, (prior_indices, prior_values) in enumerate(label_loader):
    print('Reading in label ' + str(n+1) + ' of ' + str(n_labels))
    if prior_indices.numel() == 0:
        continue
    prior_indices = torch.as_tensor(prior_indices, device=device, dtype=torch.long).squeeze()
    prior_values = torch.as_tensor(prior_values, device=device, dtype=dtype).squeeze()
    idx = tissue_index[n]
    if n == 0:
        prior = torch.sparse_coo_tensor(prior_indices[None], prior_values,
                                        [torch.Size(atlas_size).numel()]).to_dense()
        del prior_indices, prior_values
        prior = prior.reshape(torch.Size(atlas_size)).cpu().numpy()
        A[:, :, :, idx] = A[:, :, :, idx] + prior
    else:
        prior_indices = my.ind2sub(prior_indices, atlas_size)
        min_x, max_x = prior_indices[0].min().item(), prior_indices[0].max().item() + 1
        min_y, max_y = prior_indices[1].min().item(), prior_indices[1].max().item() + 1
        min_z, max_z = prior_indices[2].min().item(), prior_indices[2].max().item() + 1
        crop_atlas_size = [max_x - min_x, max_y - min_y, max_z - min_z]
        prior_indices[0] -= min_x
        prior_indices[1] -= min_y
        prior_indices[2] -= min_z
        prior = torch.sparse_coo_tensor(prior_indices, prior_values, crop_atlas_size).to_dense()
        crop = (slice(min_x, max_x), slice(min_y, max_y), slice(min_z, max_z))
        A[(*crop, idx)] = A[(*crop, idx)] + prior.cpu().numpy()
        # Hack to create maps for claustrum/reticular, LGN, and molecular layer
        if np.any(np.array([343, 368, 372, 408, 418, 562, 566, 571, 339, 354])==label_list[n]):
            MLhippo[crop] += prior.cpu().numpy()
        if np.any(np.array([102, 174])==label_list[n]):
            CLAUSTRUM[crop] += prior.cpu().numpy()
        if np.any(np.array([161, 114, 843])==label_list[n]):
            CHIASM[crop] += prior.cpu().numpy()
        if np.any(np.array([298, 305, 306, 307])==label_list[n]):
            MAMBODY[crop] += prior.cpu().numpy()
        if np.any(np.array([103, 117])==label_list[n]):
            SEPTAL[crop] += prior.cpu().numpy()
        if label_list[n]==254:
            RETICULAR[crop] += prior.cpu().numpy()
        if label_list[n]==199:
            FORNIX[crop] += prior.cpu().numpy()
        if label_list[n]==484:
            LGN[crop] += prior.cpu().numpy()
        if np.any(np.array([215,216,377,217,214,242,301,238,240,277,278,279])==label_list[n]):
            AMYGDALA[crop] += prior.cpu().numpy()

A = torch.tensor(A, dtype=dtype, device=device)
if (side=='left') or (side=='left-c') or (side=='left-ccb'):
    aff_A = np.diag([.2, .2, .2, 1])
else:
    aff_A = np.diag([-.2, .2, .2, 1])

################
# fake image
MU_CSF = 0
if (side=='left'):
    MU_WM = torch.median(Iim[Sim==2])
    MU_GM = torch.median(Iim[(Sim==3) | (Sim==17) | (Sim==18)])
    MU_WM_CEREBELLUM = torch.median(Iim[Sim==7])
    MU_GM_CEREBELLUM = torch.median(Iim[Sim==8])
    MU_CAUDATE = torch.median(Iim[Sim==11])
    MU_PUTAMEN = torch.median(Iim[Sim==12])
    MU_PALLIDUM = torch.median(Iim[Sim == 13])
else:
    MU_WM = torch.median(Iim[Sim==41])
    MU_GM = torch.median(Iim[(Sim==42) | (Sim==53) | (Sim==54)])
    MU_WM_CEREBELLUM = torch.median(Iim[Sim==46])
    MU_GM_CEREBELLUM = torch.median(Iim[Sim==47])
    MU_CAUDATE = torch.median(Iim[Sim==50])
    MU_PUTAMEN = torch.median(Iim[Sim==51])
    MU_PALLIDUM = torch.median(Iim[Sim == 52])
mid = 0.5 * MU_WM + 0.5 * MU_GM
delta = (MU_WM - MU_GM) / 16.0
MU_TH_LATERAL = mid + 2 * delta
MU_TH_MEDIAL = mid - 2 * delta
MU_RN = MU_WM + 9 * delta
MU_GM_BS = MU_WM - 1 * delta
MU_WM_BS = MU_WM + 6 * delta
MU_HYPO = mid - 3 * delta
MU_MAM_BODY = MU_WM
MU_DG_CEREBELLUM = mid - 1 * delta
MU_WM_HIPPO = mid + 1 * delta
# Kill cerebellum and brainstem if needed
if (mode=='cerebrum') or (mode=='hemi'): # We don't kill the brainstem (hard to know where to crop) and let the registration handle it
    MU_WM_CEREBELLUM = MU_GM_CEREBELLUM = MU_DG_CEREBELLUM = 0

cheating_means = torch.zeros([17], device=device, dtype=dtype)
cheating_means[0] = MU_CSF
cheating_means[1] = MU_WM
cheating_means[2] = MU_GM
cheating_means[3] = MU_WM_CEREBELLUM
cheating_means[4] = MU_GM_CEREBELLUM
cheating_means[5] = MU_CAUDATE
cheating_means[6] = MU_PUTAMEN
cheating_means[7] = MU_TH_LATERAL
cheating_means[8] = MU_TH_MEDIAL
cheating_means[9] = MU_PALLIDUM
cheating_means[10] = MU_RN
cheating_means[11] = MU_GM_BS
cheating_means[12] = MU_WM_BS
cheating_means[13] = MU_HYPO
cheating_means[14] = MU_MAM_BODY
cheating_means[15] = MU_DG_CEREBELLUM
cheating_means[16] = MU_WM_HIPPO

sigma =  torch.tensor(10.0, device=device, dtype=dtype)
if False:
    AL = torch.argmax(A, axis=-1)
    muI = cheating_means[AL]
    # sigmaI = np.sqrt(vars_ini)[AL]
    sigmaI = sigma * torch.ones(AL.shape, device=device, dtype=dtype)
else:
    muI = torch.zeros(A.shape[:-1], device=device, dtype=dtype)
    for l in range(A.shape[-1]):
        muI += (A[:,:,:,l] * cheating_means[l])
    sigmaI = sigma * torch.ones(muI.shape, device=device, dtype=dtype)

Ifake = torch.normal(muI, sigmaI)
del muI
del sigmaI
Ifake[Ifake<0] = 0
native_resolution = np.sqrt(np.sum(aff[:-1,:-1]**2, axis=0))
# we blur a bit less since the atlas is blurry already (default power factor at W/2 is 5.0)
Ifake, aff_fake = my.torch_resize(Ifake, aff_A, native_resolution, device, dtype=dtype, power_factor_at_half_width=(0.5*5.0))

# For the linear registration, we concatenate the subject-MNI transform with precomputed MNI-NextBrain transforms
if (side=='left') or (side=='left-c') or (side=='left-ccb'):
    Mmni = np.matrix('1.0959   -0.0131   -0.0233   79.9750;    0.0381    1.1359    0.0035  136.9361;  0.0735   -0.0307    1.0872   89.2956; 0 0 0 1')
else:
    Mmni = np.matrix('1.0959    0.0131    0.0233  -79.9750;   -0.0381    1.1359    0.0035  136.9361; -0.0735   -0.0307    1.0872   89.2956; 0 0 0 1')
shift_mat = np.eye(4); shift_mat[0,3]=cropping[0]; shift_mat[1,3]=cropping[1]; shift_mat[2,3]=cropping[2] # accounts for the cropping!
affine = np.linalg.inv(aff_fake) @ Mmni @ M_input_vox_to_mni_ras @ shift_mat
# Resampling
II, JJ, KK = np.meshgrid(np.arange(Iim.shape[0]), np.arange(Iim.shape[1]), np.arange(Iim.shape[2]), indexing='ij')
II = torch.tensor(II, device=device, dtype=dtype)
JJ = torch.tensor(JJ, device=device, dtype=dtype)
KK = torch.tensor(KK, device=device, dtype=dtype)
II2 = affine[0, 0] * II + affine[0, 1] * JJ + affine[0, 2] * KK + affine[0, 3]
JJ2 = affine[1, 0] * II + affine[1, 1] * JJ + affine[1, 2] * KK + affine[1, 3]
KK2 = affine[2, 0] * II + affine[2, 1] * JJ + affine[2, 2] * KK + affine[2, 3]
RSlin = my.fast_3D_interp_torch(Ifake, II2, JJ2, KK2, 'linear')
fake_filename = output_dir + '/temp_fake_linear.nii.gz'
my.MRIwrite(RSlin.detach().cpu().numpy(), aff, fake_filename)
del II, JJ, KK, II2, JJ2, KK2
del RSlin


##############################################
# OK we're ready for nonlinear registration! #
##############################################
print('Nonlinear registration of fake image with FireANTs')
# It's sad we need to write to disk... but there's no easy way of biulding a SimpleITK object
# from a generic affine vox2ras matrix :-S
ref_filename = output_dir + '/temp_reference.nii.gz'
my.MRIwrite((Iim*Mim).detach().cpu().numpy(), aff, ref_filename)

# FireANTs!
# Prepare reference image, with intensities in 1st channel, and a bunch of segmentations concatenated
image1 = Image.load_file(ref_filename, device=device_registration)
SimP = (Sim*Mim).permute([2, 1, 0])[None, None, ...]
# GM, WM (with basal forebrain), [CGM unless cerebrum mode], CA/PU/AC, TH, PA, [BG unless cerebrum mode], AM;
# Note that I had to insert hypothalamus before the amygdala (the special cases go at the end)
# In Feb'25, hypothal [13, 14], <amygdala goes here> fornix, optic chiasm, mammillary body, septal nucleus;
# I also killed WM because of ambighuities in subthalamic area (plus, its contour is modeled by cortical+subcortical GM anyway
if (side=='left') and ((mode=='invivo') or (mode=='exvivo')):
    groups = [[3,  17],#[2, 77,865],
              [8], [11, 12, 26], [10], [13], [0,  4,  5, 14, 15, 24, 41, 42, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60, 162], [819], [18], [821], [85], [843], [869]]
if (side=='left') and ((mode=='hemi') or (mode=='cerebrum')):
    groups = [[3, 17], #[2, 77,865],
              [11, 12, 26], [10], [13], [819], [18], [821], [85], [843], [869] ]
if (side=='right') and ((mode=='invivo') or (mode=='exvivo')):
    groups = [[42, 53], # [41, 77,866],
              [47], [50, 51, 58], [49], [52], [0, 43, 44, 14, 15, 24,  2,  3,  4,  5,  7,  8, 10, 11, 12, 13, 17, 18, 26, 28, 161], [820], [54], [822], [85], [844], [870]]
if (side=='right') and ((mode=='hemi') or (mode=='cerebrum')):
    groups = [[42, 53], #[41, 77,866],
              [50, 51, 58], [49], [52], [820], [54], [822], [85], [844], [870] ]
for group in groups:
    M = torch.zeros(SimP.shape, device=device_registration, dtype=dtype)
    for lab in group:
        M[SimP == lab] = 1.0
    image1.array = torch.cat([image1.array, M], dim=1)
del SimP, M

# Same for the fake image (requires resampling atlas)
image2 = Image.load_file(fake_filename, device=device_registration)
# GM, WM, [CGM unless cerebrum mode], CA/PU/AC, TH, PA, [BG unless cerebrum mode], AM (special case, see below);
# Note that I had to insert hypothalamus before the amygdala (the special cases go at the end)
# In Feb'25, hypothal [13, 14], <amygdala goes here> fornix, optic chiasm, mammillary body, septal nucleus
if (mode=='hemi') or (mode=='cerebrum'):
    frame_sets = [[2], #[1],
                  [5, 6], [7, 8], [9], [13,14]]
else:
    frame_sets = [[2], #[1],
                  [4],  [5,6], [7,8], [9], [0], [13,14]]
II, JJ, KK = np.meshgrid(np.arange(Iim.shape[0]), np.arange(Iim.shape[1]), np.arange(Iim.shape[2]), indexing='ij')
II = torch.tensor(II, device=device, dtype=dtype)
JJ = torch.tensor(JJ, device=device, dtype=dtype)
KK = torch.tensor(KK, device=device, dtype=dtype)
# fakeVox <- fakeRAS <- imageVox
# affine = np.linalg.inv(aff_A) @ my.getM(imageCOGvox[:-1,:], fakeCOGras[:-1,:]).detach().cpu().numpy()
affine = np.linalg.inv(aff_A) @ Mmni @ M_input_vox_to_mni_ras @ shift_mat
II2 = affine[0, 0] * II + affine[0, 1] * JJ + affine[0, 2] * KK + affine[0, 3]
JJ2 = affine[1, 0] * II + affine[1, 1] * JJ + affine[1, 2] * KK + affine[1, 3]
KK2 = affine[2, 0] * II + affine[2, 1] * JJ + affine[2, 2] * KK + affine[2, 3]
del II, JJ, KK


for framelist in frame_sets:
    M = torch.zeros(II2.shape, device=device, dtype=dtype)
    for frame in framelist:
        pad = 1.0 if frame==0 else 0.0
        if frame == 2:
            M += my.fast_3D_interp_torch(A[..., frame] + torch.tensor(MLhippo - AMYGDALA - SEPTAL, device=device, dtype=dtype) , II2, JJ2, KK2, 'linear', pad_value=pad)
        elif frame == 1: # should not happen anymore, but it keep it there as it does not bother me
            M += my.fast_3D_interp_torch(A[..., frame] + torch.tensor(CLAUSTRUM + RETICULAR - FORNIX - CHIASM, device=device, dtype=dtype) , II2, JJ2, KK2, 'linear', pad_value=pad)
        elif frame==6:
            M += my.fast_3D_interp_torch(A[..., frame] - torch.tensor(CLAUSTRUM, device=device, dtype=dtype), II2, JJ2, KK2, 'linear', pad_value=pad)
        elif frame==8:
            M += my.fast_3D_interp_torch(A[..., frame] - torch.tensor(LGN + RETICULAR, device=device, dtype=dtype), II2, JJ2, KK2, 'linear', pad_value=pad)
        elif frame==14:
            M += my.fast_3D_interp_torch(A[..., frame] - torch.tensor(MAMBODY, device=device, dtype=dtype), II2, JJ2, KK2, 'linear', pad_value=pad)
        else:
            M += my.fast_3D_interp_torch(A[..., frame], II2, JJ2, KK2, 'linear', pad_value=pad)
    image2.array = torch.cat([image2.array, M.permute([2, 1, 0])[None, None, ...].to(device_registration)], dim=1)
# Amygdala, fornix, chiasm, mam body, are a bit special because they are on their own
M = my.fast_3D_interp_torch(torch.tensor(AMYGDALA, device=device, dtype=dtype) , II2, JJ2, KK2, 'linear', pad_value=0)
image2.array = torch.cat([image2.array, M.permute([2, 1, 0])[None, None, ...].to(device_registration)], dim=1)
M = my.fast_3D_interp_torch(torch.tensor(FORNIX, device=device, dtype=dtype) , II2, JJ2, KK2, 'linear', pad_value=0)
image2.array = torch.cat([image2.array, M.permute([2, 1, 0])[None, None, ...].to(device_registration)], dim=1)
M = my.fast_3D_interp_torch(torch.tensor(CHIASM, device=device, dtype=dtype) , II2, JJ2, KK2, 'linear', pad_value=0)
image2.array = torch.cat([image2.array, M.permute([2, 1, 0])[None, None, ...].to(device_registration)], dim=1)
M = my.fast_3D_interp_torch(torch.tensor(MAMBODY, device=device, dtype=dtype) , II2, JJ2, KK2, 'linear', pad_value=0)
image2.array = torch.cat([image2.array, M.permute([2, 1, 0])[None, None, ...].to(device_registration)], dim=1)
M = my.fast_3D_interp_torch(torch.tensor(SEPTAL, device=device, dtype=dtype) , II2, JJ2, KK2, 'linear', pad_value=0)
image2.array = torch.cat([image2.array, M.permute([2, 1, 0])[None, None, ...].to(device_registration)], dim=1)

del II2, JJ2, KK2, M, MLhippo, CLAUSTRUM, RETICULAR, LGN, AMYGDALA, FORNIX, CHIASM, MAMBODY, SEPTAL

#FireANTs options
torch.set_grad_enabled(True)
batch1 = BatchedImages([image1])
batch2 = BatchedImages([image2])
scales = [4, 2, 1]; iterations = [300, 100, 30]
if np.mean(native_resolution)<0.7:
    scales = [6, 4, 2, 1]; iterations = [300, 100, 30, 30]
if np.mean(native_resolution)<0.5:
    scales = [8, 4, 2, 1]; iterations = [300, 100, 30, 30]
if np.mean(native_resolution)<0.3:
    scales = [12, 8, 4, 2, 1]; iterations = [300, 100, 30, 30, 30]
if np.mean(native_resolution) < 0.2:
    scales = [16, 8, 4, 2, 1]; iterations = [300, 100, 30, 30, 30]

reg = GreedyRegistration(scales=scales, iterations=iterations,
            fixed_images=batch1, moving_images=batch2,
            loss_type='custom',
            custom_loss=HybridDiceLabelDiffloss(kernel_size=args.cc_kernel_size, rel_weight_labeldiff=args.rel_weight_labeldiff),
            deformation_type='compositive',
            smooth_grad_sigma=args.smooth_grad_sigma, smooth_warp_sigma=args.smooth_warp_sigma,
            optimizer='adam', optimizer_lr=args.optimizer_lr,
            reduction='mean')  # crucial option, for proper display of loss...

start = time()
reg.optimize(save_transformed=False)
end = time()
print("Runtime nonlinear registration: ", end - start, "seconds")
moved = reg.evaluate(batch1, batch2)
reference_img = image1.itk_image
moved_image_np = moved[0, 0].detach().cpu().numpy()
moved_sitk_image = sitk.GetImageFromArray(moved_image_np)
moved_sitk_image.SetOrigin(reference_img.GetOrigin())
moved_sitk_image.SetSpacing(reference_img.GetSpacing())
moved_sitk_image.SetDirection(reference_img.GetDirection())
if args.save_atlas_nonlinear_reg:
    fake_filename_deformed = output_dir + '/atlas_nonlinear_reg.' + side + '.nii.gz'
    sitk.WriteImage(moved_sitk_image, fake_filename_deformed)
warped_coords = reg.get_warped_coordinates(batch1, batch2).to(device).detach()
del image1, image2, batch1, batch2, reg, moved, moved_sitk_image
torch.set_grad_enabled(False)

########################################################
print('Computing initial values for means and variances')
mus_ini = []
vars_ini = []
mixture_weights = []

# First, the background
if SET_BG_TO_CSF:
    x = []
    if (side=='left') or (side=='left-c') or (side=='left-ccb'):
        for l in [4, 5]: # , 24]:
            x.append(Iim_before_masking[Sim==l])
    else:
        for l in [43, 44]: # , 24]:
            x.append(Iim_before_masking[Sim==l])
    mu_bg  = torch.median(torch.concatenate(x))
    if torch.isnan(mu_bg):
        mu_bg = torch.tensor(0, dtype=dtype, device=device)
else:
    mu_bg = torch.tensor(0, dtype=dtype, device=device)

# Now, the rest of classes
for t in range(len(number_of_gmm_components)-1):
    x = []
    for l in grouping_labels[t+1]:
        x.append(Iim_before_masking[Sim==l])
    x = torch.concatenate(x) if len(x)>0 else []

    if len(x) == 0: # if there's nothing in there (eg, cerebrum mode), grab from WM or GM as needed
        nc = number_of_gmm_components[t + 1]
        if (grouping_labels[t+1][0]==7) or (grouping_labels[t+1][0]==161) or (grouping_labels[t+1][0]==46) or (grouping_labels[t+1][0]==162):
            mus_ini.append(torch.mean(mus_ini[0]) * torch.ones(nc, dtype=dtype, device=device))
            vars_ini.append(torch.mean(vars_ini[0]) * torch.ones(nc, dtype=dtype, device=device))
        elif (grouping_labels[t+1][0] == 8) or (grouping_labels[t+1][0] == 47):
            mus_ini.append(torch.mean(mus_ini[1]) * torch.ones(nc, dtype=dtype, device=device))
            vars_ini.append(torch.mean(vars_ini[1]) * torch.ones(nc, dtype=dtype, device=device))
        else:
            raise Exception('Could not initialize label group ' + str(t+1) + ' (this should not happen unless SynthSeg did really poorly)')
        mixture_weights.append((1 / float(nc)) * torch.ones(nc, dtype=dtype, device=device))
    else:

        mu = torch.median(x)
        std = 1.4826 * torch.median(torch.abs(x - mu))
        var = std ** 2
        if number_of_gmm_components[t+1]==1:
            mus_ini.append(mu[None])
            vars_ini.append(var[None])
            mixture_weights.append(torch.ones(1,dtype=dtype,device=device))
        else:
            # Estimate GMM with shared variance (avoids a component with tiny variance)
            nc = number_of_gmm_components[t+1]
            nx = len(x)
            gmm_mus = torch.linspace(mu - 0.5 * std, mu + 0.5 * std, nc, dtype=dtype, device=device)
            gmm_var= var * torch.ones(1, dtype=dtype, device=device)
            gmm_ws = (1 / float(nc)) * torch.ones(nc, dtype=dtype, device=device)
            W = torch.zeros([nx, nc], dtype=dtype, device=device)
            for its in range(200):
                # E step
                for c in range(nc):
                    W[:, c] = gmm_ws[c] / torch.sqrt(2.0 * torch.pi * torch.sqrt(gmm_var)) * torch.exp(-0.5 * (x - gmm_mus[c])**2 / gmm_var)
                normalizer = torch.sum(W + 1e-9, axis=1)
                # print(-torch.mean(torch.log(normalizer)))
                W /= normalizer[:, None]
                # M step
                denominators = torch.sum(W, axis=0)
                gmm_ws = denominators / torch.sum(denominators)
                gmm_var = 0
                for c in range(nc):
                    gmm_mus[c] = torch.sum(W[:, c] * x) / denominators[c]
                    aux = x - gmm_mus[c]
                    gmm_var += torch.sum(W[:, c] * aux * aux)
                gmm_var /= torch.sum(denominators)

            mus_ini.append(gmm_mus)
            vars_ini.append(gmm_var * torch.ones(nc, dtype=dtype, device=device))
            mixture_weights.append(gmm_ws)

mus_ini = torch.concatenate(mus_ini)
vars_ini = torch.concatenate(vars_ini)
mixture_weights = torch.concatenate(mixture_weights)
# Replace -1 variance with minimum otherwise
vars_ini[vars_ini < 0] = torch.min(vars_ini[vars_ini>0])

########################################################

# EM for GMM parameters
print('Estimating GMM parameters')

print('  Resizing image, mask, and coordinates')
I_r, aff_r = my.torch_resize(Iim_before_masking, aff, resolution, device, dtype=dtype)
M_r, _ = my.torch_resize(Mim, aff, resolution, device, dtype=dtype)
# smoothen M_r a bit
kernel = torch.zeros([3,3,3], device=device, dtype=dtype)
kernel[1,1,:] = 0.1; kernel[1,:,1] = 0.1; kernel[:,1,1] = 0.1; kernel[1,1,1] = 0.4
for _ in range(3):
    M_r = torch.conv3d(M_r[None, None, ...], kernel[None, None, ...], bias=None, stride=1, padding=[1, 1, 1]).squeeze()
I_r[M_r<0.5] = mu_bg
Iim_shape = Iim.shape
del Iim, Iim_before_masking, Mim, Sim

print('  Deforming atlas')
# OK so we need to create a field of coordinates and then resize it.
# The tricky bit is concatenating the nonlinear transform (in [0,1]) with the linear
wc_r, _ = my.torch_resize(warped_coords[0].permute([2,1,0,3]), aff, resolution, device, dtype=dtype)
del warped_coords
# we'll worry about skip later
# first bit: from [-1,1] of linearly registered, to absolute coordinates
T1 = torch.eye(4, device=device, dtype=dtype)
for k in range(3):
    T1[k, k] = 0.5 * (Iim_shape[k] - 1)
    T1[k, -1] = T1[k, k]
# second bit: vox2vox transform to fake  space  fake_vox <- image_vox
# T2 = my.getM(imageCOGvox[:-1,:], fakeCOGvox[:-1,:])
T2 = torch.tensor(np.linalg.inv(aff_fake) @ Mmni @ M_input_vox_to_mni_ras @ shift_mat, device=device, dtype=dtype) # TODO: is this correct?
# third bit: from vox to [-1, 1] coordinates that can be used with atlas at any resolution
T3 = torch.eye(4, device=device, dtype=dtype)
for k in range(3):
    T3[k, k] = 2 / (Ifake.shape[k] - 1)
    T3[k, -1] = -1
T = T3 @ T2 @ T1
I = T[0, 0] * wc_r[..., 0] + T[0, 1] * wc_r[..., 1] + T[0, 2] * wc_r[..., 2] + T[0, 3]
J = T[1, 0] * wc_r[..., 0] + T[1, 1] * wc_r[..., 1] + T[1, 2] * wc_r[..., 2] + T[1, 3]
K = T[2, 0] * wc_r[..., 0] + T[2, 1] * wc_r[..., 1] + T[2, 2] * wc_r[..., 2] + T[2, 3]
del wc_r
# we can now resample
priors = grid_sample(A.permute([3,0,1,2])[None, ...],
                     torch.stack([K[::skip,::skip,::skip],
                                  J[::skip,::skip,::skip],
                                  I[::skip,::skip,::skip]], axis=-1)[None,...], align_corners=True)
priors = torch.permute(priors[0], [1, 2, 3, 0])
# Deal with voxels outside the FOV
missing_mass = 1 - torch.sum(priors, axis=-1)
priors[..., 0] += missing_mass

####
print('  EM for parameter estimation')
# data
means = torch.tensor([mu_bg, *mus_ini], device=device, dtype=dtype)
var_bg = torch.min(vars_ini)
variances = torch.tensor([var_bg, *vars_ini], device=device, dtype=dtype)
weights = torch.tensor([1.0, *mixture_weights], dtype=dtype, device=device)
W = torch.zeros([*priors.shape[:-1], number_of_gmm_components.sum()], dtype=dtype, device=device)
x = I_r[::skip, ::skip, ::skip]

# We now put a Scaled inverse chi-squared prior on the variances to prevent them from going to zero
prior_count = 100 / (resolution ** 3) / (skip ** 3)
prior_variance = var_bg
loglhood_old = -10000
eps = 1e-7
for em_it in range(100):
    # E step
    for c in range(n_tissues):
        prior = priors[:, :, :, c]
        num_components = number_of_gmm_components[c]
        for g in range(num_components):
            gaussian_number = sum(number_of_gmm_components[:c]) + g
            d = x - means[gaussian_number]
            W[:, :, :, gaussian_number] = weights[gaussian_number] * prior * torch.exp(
                -d * d / (2 * variances[gaussian_number])) / torch.sqrt(
                2.0 * torch.pi * variances[gaussian_number])

    normalizer = eps + torch.sum(W, dim=-1, keepdim=True)
    loglhood = torch.mean(torch.log(normalizer)).detach().cpu().numpy()
    W = W / normalizer

    # M step
    prior_loglhood = torch.zeros(1,dtype=dtype,device=device)
    for c in range(number_of_gmm_components.sum()):
        # crucially, we skip the background when we update the parameters (but we still add it to the cost)
        if c > 0:
            norm = eps + torch.sum(W[:, :, :, c])
            means[c] = torch.sum(x * W[:, :, :, c]) / norm
            d = x - means[c]
            variances[c] = (torch.sum(d * d * W[:, :, :, c]) + prior_count * prior_variance) / ( norm + prior_count + 2)
        v = variances[c]
        prior_loglhood = prior_loglhood - ((1 + 0.5 * prior_count) * torch.log(v) + 0.5 * prior_count * prior_variance / v) / torch.numel(normalizer)
    loglhood = loglhood + prior_loglhood.detach().cpu().numpy()

    mixture_weights = torch.sum(W[:, :, :, 1:].reshape([np.prod(priors.shape[:-1]), number_of_gmm_components.sum() - 1]) + eps, axis=0)
    for c in range(n_tissues - 1):
        # mixture weights are normalized (those belonging to one mixture sum to one)
        num_components = number_of_gmm_components[c + 1]
        gaussian_numbers = torch.tensor(np.sum(number_of_gmm_components[1:c + 1]) + \
                                        np.array(range(num_components)), device=device, dtype=dtype).long()

        mixture_weights[gaussian_numbers] /= torch.sum(mixture_weights[gaussian_numbers])

    weights[1:] = mixture_weights

    if (torch.sum(torch.isnan(means)) > 0) or (torch.sum(torch.isnan(variances)) > 0):
        print('nan in Gaussian parameters...')
        import pdb;

        pdb.set_trace()

    print('         Step %d of EM, -loglhood = %.6f' % (em_it + 1, -loglhood ), flush=True)
    if (loglhood - loglhood_old) < TOL:
        print('         Decrease in loss below tolerance limit')
        break
    else:
        loglhood_old = loglhood


###########

print('Computing Gaussians at full resolution (we will reuse over and over)')
GAUSSIAN_LHOODS = torch.zeros([*I_r.shape, sum(number_of_gmm_components)], dtype=dtype, device=device)
for c in range(sum(number_of_gmm_components)):
    # The 1e-9 ensures no zeros were prior is not zero
    GAUSSIAN_LHOODS[..., c] = 1.0 / torch.sqrt(2 * math.pi * variances[c]) * torch.exp(
    -0.5 * torch.pow(I_r - means[c], 2.0) / variances[c]) + eps

print('Computing normalizers (faster to do now with clustered priors)')
# We deform one class at the time; slower, but less memory
normalizers = torch.zeros(GAUSSIAN_LHOODS.shape[:-1], dtype=dtype, device=device)
# normalizers[h] = eps
gaussian_number = 0
for c in range(A.shape[-1]):
    prior = grid_sample(A[None, None, ..., c], torch.stack([K, J, I], axis=-1)[None, ...], align_corners=True)[0,0,...]
    if c==0: # background
        prior[(I < (-1)) | (I > 1) | (J < (-1)) | (J > 1) | (K < (-1)) | (K > 1)] = 1.0
    lhood = torch.zeros_like(prior)
    for g in range(number_of_gmm_components[c]):
        lhood += (weights[gaussian_number] * GAUSSIAN_LHOODS[..., gaussian_number])
        gaussian_number += 1
    normalizers += (prior * lhood)
Ashape = A.shape
del A

print('Deforming one label at the time')
names, colors = my.read_LUT(LUT_file)
seg = torch.zeros(normalizers.shape, dtype=torch.int, device=device)
seg_rgb = torch.zeros([*normalizers.shape, 3], dtype=dtype, device=device)
max_p = torch.zeros(normalizers.shape, dtype=dtype, device=device)
vols = torch.zeros(n_labels, device=device, dtype=dtype)

# TODO: choose good number of workers/prefetch factor
for n, (prior_indices, prior_values) in enumerate(label_loader):
    print('Deforming label ' + str(n + 1) + ' of ' + str(n_labels))

    if prior_indices.numel() == 0:
        continue
    prior_indices = torch.as_tensor(prior_indices, device=device, dtype=torch.long).squeeze()
    prior_values = torch.as_tensor(prior_values, device=device, dtype=dtype).squeeze()

    if n == 0:
        # background
        prior = torch.sparse_coo_tensor(prior_indices[None], prior_values,
                                        [torch.Size(atlas_size).numel()]).to_dense()
        del prior_indices, prior_values
        prior = prior.reshape(torch.Size(atlas_size))

    else:

        # find bounding box of label in atlas space
        prior_indices = my.ind2sub(prior_indices, atlas_size)
        min_x, max_x = prior_indices[0].min().item(), prior_indices[0].max().item()
        min_y, max_y = prior_indices[1].min().item(), prior_indices[1].max().item()
        min_z, max_z = prior_indices[2].min().item(), prior_indices[2].max().item()
        crop_atlas_size = [max_x - min_x + 1, max_y - min_y + 1, max_z - min_z + 1]
        prior_indices[0] -= min_x
        prior_indices[1] -= min_y
        prior_indices[2] -= min_z
        prior = torch.sparse_coo_tensor(prior_indices, prior_values, crop_atlas_size).to_dense()
        del prior_indices, prior_values

    skip_this_label = False
    if n==0:
        Irescaled = I
        Jrescaled = J
        Krescaled = K
        lr_crop = (slice(None),) * 3
    else:
        # find bounding box of label in MRI space
        Irescaled = I * ((Ashape[0] - 1) / (max_x - min_x)) + ( (Ashape[0] - 1 - min_x - max_x) / (max_x - min_x) )
        Jrescaled = J * ((Ashape[1] - 1) / (max_y - min_y)) + ( (Ashape[1] - 1 - min_y - max_y) / (max_y - min_y) )
        Krescaled = K * ((Ashape[2] - 1) / (max_z - min_z)) + ( (Ashape[2] - 1 - min_z - max_z) / (max_z - min_z) )
        mask = (Irescaled >= (-1))
        mask &= (Irescaled <= 1)
        mask &= (Jrescaled >= (-1))
        mask &= (Jrescaled <= 1)
        mask &= (Krescaled >= (-1))
        mask &= (Krescaled <= 1)
        if mask.any()==False:
            skip_this_label = True
        else:
            nx, ny, nz = mask.shape
            tmp = mask.reshape([nx, -1]).any(-1).nonzero()
            lr_min_x, lr_max_x = tmp.min().item(), tmp.max().item() + 1
            tmp = mask.movedim(0, -1).reshape([ny, -1]).any(-1).nonzero()
            lr_min_y, lr_max_y = tmp.min().item(), tmp.max().item() + 1
            tmp = mask.reshape([-1, nz]).any(0).nonzero()
            lr_min_z, lr_max_z = tmp.min().item(), tmp.max().item() + 1
            del tmp, mask
            lr_crop = (slice(lr_min_x, lr_max_x), slice(lr_min_y, lr_max_y), slice(lr_min_z, lr_max_z))

    if skip_this_label==False:

        prior_resampled = grid_sample(prior[None, None, ...], torch.stack([Krescaled[lr_crop], Jrescaled[lr_crop], Irescaled[lr_crop]], axis=-1)[None, ...], align_corners=True)[0, 0, ...]
        if n==0: # background
            prior_resampled[(Krescaled[lr_crop]<(-1)) | (Krescaled[lr_crop]>1) | (Jrescaled[lr_crop]<(-1)) | (Jrescaled[lr_crop]>1)  | (Irescaled[lr_crop]<(-1)) | (Irescaled[lr_crop]>1)  ] = 1.0
        del Irescaled; del Jrescaled; del Krescaled
        num_components = number_of_gmm_components[tissue_index[n]]
        gaussian_numbers = torch.tensor(np.sum(number_of_gmm_components[:tissue_index[n]]) + \
                                np.array(range(num_components)), device=device, dtype=dtype).long()
        lhood = torch.sum(GAUSSIAN_LHOODS[:, :, :, gaussian_numbers] * weights[None, None, None, gaussian_numbers], 3)
        post = torch.squeeze(prior_resampled)
        post *= lhood[lr_crop]
        post /= normalizers[lr_crop]
        if n==0:
            post[torch.isnan(post)] = 1.0
        else:
            post[torch.isnan(post)] = 0.0
            if (mode=='hemi') or (mode=='cerebrum'): # cerebrum only: kill labels outside mask
                post *= M_r[lr_crop]
        # Also, for cerebrum only, kill cerebellar labels (not brainstem as we don't know where to stop)
        if ((mode == 'hemi') or (mode == 'cerebrum ')) and (np.sum(np.array([595, 597, 715, 721, 751, 752, 846])==label_list[n])):
            post[:] = 0
        del prior_resampled
        vols[n] = torch.sum(post) * (resolution ** 3)
        mask = (post > max_p[lr_crop])
        max_p[lr_crop][mask] = post[mask]
        lab = int(label_list[n])
        seg[lr_crop].masked_fill_(mask, lab)
        del mask
        for c in range(3):
            seg_rgb[(*lr_crop, c)].add_(post, alpha=colors[lab][c])
print('\n')

########################################################

print('Writing results to disk')
my.MRIwrite(seg.detach().cpu().numpy().astype(np.uint16), aff_r, output_dir + '/seg.' + side + '.nii.gz')
if args.write_rgb:
    my.MRIwrite(seg_rgb.detach().cpu().numpy().astype(np.uint8), aff_r, output_dir + '/seg.' + side + '.rgb.nii.gz')
vols = vols.detach().cpu().numpy()
with open(output_dir + '/vols.' + side + '.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    aux = label_list[1:]
    row = []
    for l in aux:
        row.append(names[int(l)])
    writer.writerow(row)
    row = []
    for j in range(1, len(vols)):
        row.append(str(vols[j]))
    writer.writerow(row)
# Copy LUT file, for convenience
a = os.system('cp ' + LUT_file + ' ' + output_dir +  '/lut.txt >/dev/null')
if a==0:
    LUT_file = output_dir +  '/lut.txt'

# Clean up
os.system('rm -rf '  + output_dir +  '/temp*.nii.gz >/dev/null')

# Print commands to visualize output
print('You can try these commands to visualize outputs:')
cmd = '  freeview -v ' + input_volume + ' -v ' + output_dir + '/supersynth.nii.gz:colormap=lut '
if args.write_bias_corrected and (skip_bf==False):
    cmd = cmd + ' -v ' + output_dir + '/bias.corrected.nii.gz '
cmd = cmd + ' -v ' + output_dir + '/seg.' + side + '.nii.gz:colormap=lut:lut=' + LUT_file
if args.write_rgb:
    cmd = cmd + ' -v ' + output_dir + '/seg.' + side + '.rgb.nii.gz:rgb=true'
if args.save_atlas_nonlinear_reg:
    cmd = cmd + ' ' + fake_filename_deformed
print(cmd)
print('  oocalc ' + output_dir + '/vols.' + side + '.csv')
print(' ')
print('All done!')

##################

now2 = datetime.now()

current_time = now2.strftime("%H:%M:%S")
print("Current Time =", current_time)

runtime = now2 - now

print("Running Time =", runtime)
