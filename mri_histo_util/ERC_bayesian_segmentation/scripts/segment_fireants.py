#Imports
import os
import sys
basepath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(basepath)
from pathlib import Path
BASE_PATH = str(Path(__file__).resolve().parents[1])
sys.path.append(str(Path(__file__).resolve().parents[1]))
FIREANTS_PATH = os.path.join(os.environ.get('FREESURFER_HOME'),'python/packages/ERC_bayesian_segmentation/ext')
sys.path.append(FIREANTS_PATH)
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
#from fireants.io import Image,BatchedImages
#from fireants.registration import AffineRegistration, GreedyRegistration
#from ext.fireants.io import Image,BatchedImages
#from ext.fireants.registration import AffineRegistration, GreedyRegistration
from ext.fireants_trimmed import Image, BatchedImages
from ext.fireants_trimmed import GreedyRegistration
import SimpleITK as sitk
from time import time
#from ERC_bayesian_segmentation.custom_cc_labeldiff_loss import HybridDiceLabelDiffloss
from ext.fireants_trimmed import HybridDiceLabelDiffloss


########################################################

parser = argparse.ArgumentParser(description='Bayesian segmentation.')
parser.add_argument("--i", help="Image to segment.")
parser.add_argument("--i_seg", help="SynthSeg of image to segment (must include parcels). Will be computed if it does not exist")
parser.add_argument("--i_field", help="Registration to MNI, provided by EasyReg, for left-right division etc. Will be computed if if does not exist")
parser.add_argument("--atlas_dir", help="Atlas directory")
parser.add_argument("--side", help="Hemisphere to segment (left or right).")
parser.add_argument("--bf_mode", help="bias field basis function: dct, polynomial, or hybrid", default="dct")
parser.add_argument("--o", help="Output directory.")
parser.add_argument("--write_rgb", action="store_true", help="Write soft segmentation to dis as RGB file.")
parser.add_argument("--write_bias_corrected", action="store_true", help="Write bias field corrected image to disk")
parser.add_argument("--cpu", action="store_true", help="Use CPU instead of GPU")
parser.add_argument("--threads", type=int, default=1, help="(optional) Number of CPU cores to be used. Default is 1. You can use -1 to use all available cores")
parser.add_argument("--skip", type=int, default=1, help="(optional) Skipping factor to easy memory requirements of priors when estimating Gaussian parameters. Default is 1.")
parser.add_argument("--resolution", type=float, default=0.4, help="(optional) Resolution of output segmentation")
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
if args.i_seg is None:
    raise Exception('SynthSeg file  of input image is required')
if args.i_field is None:
    raise Exception('Atlas registration file of input image is required')
if args.atlas_dir is None:
    raise Exception('Atlas directory must be provided')
if args.o is None:
    raise Exception('Output directory must be provided')
if args.side is None:
    raise Exception('side must be provided (left or right)')
else:
    if (args.side!='left') and (args.side!='right'):
        raise Exception('Side must be left or right but you specified: ' + args.side)

########################################################

if args.cpu:
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
else:
    os.environ["CUDA_VISIBLE_DEVICES"] = "0"

########################################################

# limit the number of threads to be used if running on CPU
if args.threads<0:
    args.threads = os.cpu_count()
    print('using all available threads ( %s )' % args.threads)
else:
    print('using %s thread(s)' % args.threads)
torch.set_num_threads(args.threads)

########################################################

# Disable gradients
torch.set_grad_enabled(False)


################

# Input data
input_volume = args.i
input_seg = args.i_seg
atlas_field = args.i_field
atlas_dir = args.atlas_dir
LUT_file = os.path.join(BASE_PATH, 'data_simplified', 'AllenAtlasLUT')
output_dir = args.o
skip_bf = args.skip_bf
bf_mode = args.bf_mode
side = args.side
resolution = args.resolution
if (resolution<=0):
    print('Resolution must be non-negative; Exitting...')
    sys.exit(1)
skip = args.skip
if (skip<1):
    print('Skip cannot be less than 1; exitting...')
    sys.exit(1)

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

if torch.cuda.is_available():
    print('Using the GPU')
    device = torch.device('cuda:0')
else:
    print('Using the CPU')
    device = torch.device('cpu')

########################################################

now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

########################################################
print('Reading input image')
Iim, aff = my.MRIread(input_volume)
Iim = np.squeeze(Iim)

########################################################
synthseg_exists = False
if os.path.isfile(input_seg):
    print('Found input synthseg segmentation')
    synthseg_exists = True

if synthseg_exists:
    print('Input segmentation exists; making sure it includes parcellation!')
    tmp, _ = my.MRIread(input_seg)
    if np.sum(tmp>1000)==0:
        raise Exception('Segmentation does not include parcellation! Please use different file or re-run SynthSeg with --parc')
else:
    print('Running SynthSeg')
    cmd = 'mri_synthseg --i ' + input_volume +  ' --o '  + input_seg + ' --threads ' + str(args.threads) + ' --cpu --parc --robust '
    cmd = cmd + ' --vol ' + output_dir + '/SynthSeg_volumes.csv'
    a = os.system(cmd + ' >/dev/null')
    if a > 0:
        print('Error in mri_synthseg; exitting...')
        sys.exit(1)


TMP_RESAMPLED = output_dir + '/tmp.resampled.mgz'
a = os.system('mri_convert ' + input_seg + ' ' + TMP_RESAMPLED + ' -rl ' + input_volume + ' -rt nearest -odt float >/dev/null')
if a>0:
    print('Error in mri_convert; exitting...')
    sys.exit(1)
Sim, _ = my.MRIread(TMP_RESAMPLED)
os.system('rm -rf ' + TMP_RESAMPLED)

########################################################
if skip_bf==False:
    print('Correcting bias field')
    print('   Trying model with polynomial basis functions')
    try:
        Iim, _ = bf.correct_bias(Iim, Sim, maxit=100, penalty=0.1, order=6, device=device, dtype=dtype, basis=bf_mode)
    except:
        if device.type=='cpu':
            raise Exception('Bias correction failed (out of memory?)')
        else:
            print('Bias correction on GPU failed; trying with CPU')
            Iim, _ = bf.correct_bias(Iim, Sim, maxit=100, penalty=0.1, order=4, device='cpu', dtype=dtype, basis=bf_mode)
    if args.write_bias_corrected:
        aux = (Iim / np.max(Iim) * 255).astype(np.uint8)
        my.MRIwrite(aux, aff, output_dir + '/bias.corrected.nii.gz')

print('Normalizing intensities')
Iim = Iim * 110 / np.median(Iim[(Sim==2) | (Sim==41)])

# We should do tensors at this point...
Sim = torch.tensor(Sim, dtype=torch.int, device=device)
Iim = torch.tensor(Iim, dtype=dtype, device=device)

########################################################

# TODO: do this within python at some point...
print('Subdividing brainstem into left and right halves, and cropping bottom if needed')
BRAINSTEM = 16
BRAINSTEM_L = 161
BRAINSTEM_R = 162

if os.path.isfile(atlas_field):
    print('   Registration file found; no need to run EasyReg!')
else:
    print('   Running EasyReg')
    mni_image = os.path.join(BASE_PATH, 'data_mni', 'mni.nii.gz')
    mni_seg = os.path.join(BASE_PATH, 'data_mni', 'mni.synthseg.nii.gz')
    a = os.system('mri_easyreg --flo ' + mni_image +
                  ' --ref ' + input_volume + ' --flo_seg ' + mni_seg + ' --ref_seg ' + input_seg +
                  ' --threads ' + str(args.threads) + ' --fwd_field ' + atlas_field + '  >/dev/null')
    if a>0:
        print('Error in mri_easyreg; exitting...')
        sys.exit(1)

FIELD, FIELDaff = my.MRIread(atlas_field)
FIELD = torch.tensor(FIELD, dtype=dtype, device=device)
# sanity check
if np.sum(np.abs(aff-FIELDaff))>1e-6:
    raise Exception('affine matrix of deformation field does not coincide with that of input')
    sys.exit(1)

# kill bottom of medulla
CSF = 24
Sim[((Sim==BRAINSTEM) | (Sim==CSF)) & (FIELD[...,2]<(-60))] = 0
# Subdivide
LEFT = (FIELD[...,0]<0)
Sim[(Sim==BRAINSTEM) & LEFT] = BRAINSTEM_L
Sim[(Sim==BRAINSTEM) & (LEFT==0)] = BRAINSTEM_R
del LEFT, FIELD

#######################################

# Prepare data for hemisphere at hand
print('  Creating and applying mask for ' + side + ' hemisphere')
if side=='left':
    M = ( (Sim < 30) | (Sim==161) | ( (Sim > 1000) & (Sim < 2000) ) )
else: # right hemi
    M = ( ( (Sim > 40) & (Sim<100) )  | (Sim==162) | (Sim > 2000) )
if True:
    M[Sim==4] = 0
    M[Sim==5] = 0
    M[Sim==43] = 0
    M[Sim==44] = 0
M[Sim==14] = 0
M[Sim==15] = 0
M[Sim==24] = 0
M[Sim==0] = 0
# I now do this with the resampled mask later on, to avoid blurring mask edges
# Iim[M==0] = 0
Mim, cropping = my.cropLabelVolTorch(M, margin=5)
del M
Sim = my.applyCropping(Sim, cropping)
Iim = my.applyCropping(Iim, cropping)
aff[:3, -1] = aff[:3, -1] + aff[:-1, :-1] @ np.array([cropping[0].detach().cpu().numpy(), cropping[1].detach().cpu().numpy(), cropping[2].detach().cpu().numpy()])

##### Version with and without parcels
SimParc = Sim.clone()
Sim[Sim >= 2000] = 42
Sim[Sim > 1000] = 3


########################################################

# Read atlas
print('Reading in atlas')

# ASEG labels are hard coded for each hemi, no biggie
if side=='left':
    aseg_label_list = np.array(
        [2, 7, 8, 10, 11, 12, 13, 17, 18, 26, 28, 161, 1001, 1002, 1003, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012,
         1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030,
         1031, 1032, 1033, 1034, 1035]).astype(int)
else:
    aseg_label_list = np.array(
        [41, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60, 162, 2001, 2002, 2003, 2005, 2006, 2007, 2008, 2009, 2010, 2011,
         2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026, 2027, 2028, 2029,
         2030, 2031, 2032, 2033, 2034, 2035]).astype(int)

# Get the label groupings and atlas labels from the config files
tissue_index, grouping_labels, label_list, number_of_gmm_components = relab.get_tissue_settings(
            os.path.join(BASE_PATH, 'data_simplified', 'atlas_names_and_labels.yaml'),
            os.path.join(BASE_PATH, 'data_simplified', 'combined_atlas_labels_fireants.yaml'),
            os.path.join(BASE_PATH, 'data_simplified', 'combined_aseg_labels_fireants.yaml'),
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
prefetch = 4
workers = 2
prefetch_factor = max(prefetch//workers, 1)
label_loader = DataLoader(LabelDataset(atlas_names), num_workers=workers, prefetch_factor=prefetch_factor)
A = np.zeros([*atlas_size, n_tissues], dtype=numpy_dtype)
# We keep track of these probability masses we use to correct differences in labeling between A and SynthSeg
MLhippo = np.zeros(atlas_size, dtype=numpy_dtype)
CLAUSTRUM = np.zeros(atlas_size, dtype=numpy_dtype)
RETICULAR = np.zeros(atlas_size, dtype=numpy_dtype)
LGN =  np.zeros(atlas_size, dtype=numpy_dtype)
AMYGDALA =  np.zeros(atlas_size, dtype=numpy_dtype)
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
        if label_list[n]==254:
            RETICULAR[crop] += prior.cpu().numpy()
        if label_list[n]==484:
            LGN[crop] += prior.cpu().numpy()
        if np.any(np.array([215,216,377,217,214,242,301,238,240,277,278,279])==label_list[n]):
            AMYGDALA[crop] += prior.cpu().numpy()

A = torch.tensor(A, dtype=dtype, device=device)
if side=='left':
    aff_A = np.diag([.2, .2, .2, 1])
else:
    aff_A = np.diag([-.2, .2, .2, 1])

################
# fake image
MU_CSF = 0
if side=='left':
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

# For the linear registration, we used the centroids, as in EasyReg.
if side=='left':
    labels = np.array([2,4,5,7,8,10,11,12,13,161,17,18,26,28,1001,1002,1003,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1034,1035])
else:
    labels = np.array([41,43,44,46,47,49,50,51,52,162, 53, 54, 58, 60, 2001,2002,2003,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025,2026,2027,2028,2029,2030,2031,2032,2033,2034,2035])
nlab = len(labels)
imageCOGvox = torch.zeros([4, nlab], device=device, dtype=dtype)
ok = torch.ones(nlab, device=device, dtype=torch.bool)
for l in range(nlab):
    aux = torch.where(SimParc == labels[l])
    if len(aux[0]) > 50:
        imageCOGvox[0, l] = torch.median(aux[0].to(dtype))
        imageCOGvox[1, l] = torch.median(aux[1].to(dtype))
        imageCOGvox[2, l] = torch.median(aux[2].to(dtype))
        imageCOGvox[3, l] = 1.0
    else:
        ok[l] = False
imageCOGras = torch.tensor(aff, device=device, dtype=dtype) @ imageCOGvox
fakeCOGras = torch.tensor([[49.4,63.4,44.4,59.4,51.4,67.4,63.4,51.4,58.4,73.4,52.4,56.4,71.4,70.4,18.4,76.4,38.4,70.4,62.4,42.4,32.4,24.4,74.4,46.4,55.4,61.4,75.4,15.4,54.4,75.4,27.4,33.4,27.4,66.4,28.4,77.4,31.4,73.4,75.4,42.4,70.4,58.4,22.4,18.4,71.4,53.4,28.4,39.4],[117.4,107.4,116.4,81.4,69.4,113.4,147.4,137.4,131.4,105.4,112.4,130.4,147.4,117.4,86.4,161.4,147.4,41.4,131.4,81.4,59.4,98.4,79.4,31.4,169.4,53.4,179.4,111.4,99.4,103.4,155.4,181.4,171.4,39.4,112.4,113.4,127.4,63.4,183.4,189.4,168.4,62.4,128.4,95.4,206.4,152.4,111.4,136.4],[110.4,103.4,70.4,46.4,50.4,93.4,99.4,85.4,84.4,53.4,67.4,63.4,77.4,75.4,95.4,122.4,144.4,111.4,50.4,67.4,124.4,55.4,102.4,95.4,64.4,84.4,67.4,68.4,67.4,148.4,98.4,69.4,89.4,102.4,138.4,129.4,138.4,128.4,87.4,109.4,143.4,147.4,77.4,122.4,68.4,42.4,94.4,83.4],[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]], device=device, dtype=dtype)
if side=='right':
    fakeCOGras[0, :] = -fakeCOGras[0, :]
fakeCOGvox = torch.linalg.inv(torch.tensor(aff_fake, device=device, dtype=dtype)) @ fakeCOGras
imageCOGvox = imageCOGvox[:,ok]
imageCOGras = imageCOGras[:,ok]
fakeCOGras = fakeCOGras[:,ok]
fakeCOGvox = fakeCOGvox[:,ok]
# Resampling
II, JJ, KK = np.meshgrid(np.arange(Iim.shape[0]), np.arange(Iim.shape[1]), np.arange(Iim.shape[2]), indexing='ij')
II = torch.tensor(II, device=device, dtype=dtype)
JJ = torch.tensor(JJ, device=device, dtype=dtype)
KK = torch.tensor(KK, device=device, dtype=dtype)
affine = my.getM(imageCOGvox[:-1,:], fakeCOGvox[:-1,:])
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
image1 = Image.load_file(ref_filename, device=device)
SimP = (Sim*Mim).permute([2, 1, 0])[None, None, ...]
# GM, WM, CGM, CA/PU/AC, TH, PA, BG, AM
if side=='left':
    groups = [[3,  17],  [2],  [8], [11, 12, 26], [10], [13], [0,  4,  5, 14, 15, 24, 41, 42, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60, 162], [18]]
else:
    groups = [[42, 53], [41], [47], [50, 51, 58], [49], [52], [0, 43, 44, 14, 15, 24,  2,  3,  4,  5,  7,  8, 10, 11, 12, 13, 17, 18, 26, 28, 161], [54]]
for group in groups:
    M = torch.zeros(SimP.shape, device=device, dtype=dtype)
    for lab in group:
        M[SimP == lab] = 1.0
    image1.array = torch.cat([image1.array, M], dim=1)
del SimP, M

# Same for the fake image (requires resampling atlas)
image2 = Image.load_file(fake_filename, device=device)
# GM, WM, CGM, CA/PU/AC, TH, PA, BG, AM
frame_sets = [[2], [1], [4],  [5,6], [7,8], [9], [0]]
II, JJ, KK = np.meshgrid(np.arange(Iim.shape[0]), np.arange(Iim.shape[1]), np.arange(Iim.shape[2]), indexing='ij')
II = torch.tensor(II, device=device, dtype=dtype)
JJ = torch.tensor(JJ, device=device, dtype=dtype)
KK = torch.tensor(KK, device=device, dtype=dtype)
affine = np.linalg.inv(aff_A) @ my.getM(imageCOGvox[:-1,:], fakeCOGras[:-1,:]).detach().cpu().numpy()
II2 = affine[0, 0] * II + affine[0, 1] * JJ + affine[0, 2] * KK + affine[0, 3]
JJ2 = affine[1, 0] * II + affine[1, 1] * JJ + affine[1, 2] * KK + affine[1, 3]
KK2 = affine[2, 0] * II + affine[2, 1] * JJ + affine[2, 2] * KK + affine[2, 3]
del II, JJ, KK

for framelist in frame_sets:
    M = torch.zeros(II2.shape, device=device, dtype=dtype)
    for frame in framelist:
        pad = 1.0 if frame==0 else 0.0
        if frame == 2:
            M += my.fast_3D_interp_torch(A[..., frame] + torch.tensor(MLhippo - AMYGDALA, device=device, dtype=dtype) , II2, JJ2, KK2, 'linear', pad_value=pad)
        elif frame == 1:
            M += my.fast_3D_interp_torch(A[..., frame] + torch.tensor(CLAUSTRUM + RETICULAR, device=device, dtype=dtype) , II2, JJ2, KK2, 'linear', pad_value=pad)
        elif frame==6:
            M += my.fast_3D_interp_torch(A[..., frame] - torch.tensor(CLAUSTRUM, device=device, dtype=dtype), II2, JJ2, KK2, 'linear', pad_value=pad)
        elif frame==8:
            M += my.fast_3D_interp_torch(A[..., frame] - torch.tensor(LGN + RETICULAR, device=device, dtype=dtype), II2, JJ2, KK2, 'linear', pad_value=pad)
        else:
            M += my.fast_3D_interp_torch(A[..., frame], II2, JJ2, KK2, 'linear', pad_value=pad)
    image2.array = torch.cat([image2.array, M.permute([2, 1, 0])[None, None, ...]], dim=1)
# Amygdala is a bit special because it's on its own
M = my.fast_3D_interp_torch(torch.tensor(AMYGDALA, device=device, dtype=dtype) , II2, JJ2, KK2, 'linear', pad_value=0)
image2.array = torch.cat([image2.array, M.permute([2, 1, 0])[None, None, ...]], dim=1)
del II2, JJ2, KK2, M, MLhippo, CLAUSTRUM, RETICULAR, LGN, AMYGDALA

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
warped_coords = reg.get_warped_coordinates(batch1, batch2).detach()
del image1, image2, batch1, batch2, reg, moved, moved_sitk_image
torch.set_grad_enabled(False)
# if False: # this is how FireANTs does it
#     resampled = grid_sample(RSlin.permute([2,1,0])[None, None, ...], warped_coords, mode='bilinear', align_corners=True).squeeze().permute([2,1,0])
#     my.MRIwrite(torch.squeeze(resampled).detach().cpu().numpy(),aff,'/tmp/test.nii.gz')
# if False: # an easier way to look at it?
#     resampled = grid_sample(RSlin[None, None, ...], warped_coords.flip(4), mode='bilinear', align_corners=True).squeeze().permute([2,1,0])
#     my.MRIwrite(torch.squeeze(resampled).detach().cpu().numpy(),aff,'/tmp/test2.nii.gz')
# if False: # yet another way to look at it?
#     resampled = grid_sample(RSlin[None, None, ...], warped_coords.permute([0,3,2,1,4]).flip(4), mode='bilinear', align_corners=True).squeeze()
#     my.MRIwrite(torch.squeeze(resampled).detach().cpu().numpy(),aff,'/tmp/test3.nii.gz')


########################################################
print('Computing initial values for means and variances')
mus_ini = []
vars_ini = []
mixture_weights = []

for t in range(len(number_of_gmm_components)-1):
    x = []
    for l in grouping_labels[t+1]:
        x.append(Iim[SimParc==l])

    if len(x) > 0:
        x = torch.concatenate(x)
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

if SET_BG_TO_CSF:
    x = []
    if side=='left':
        for l in [4, 5]: # , 24]:
            x.append(Iim[Sim==l])
    else:
        for l in [43, 44]: # , 24]:
            x.append(Iim[Sim==l])
    mu_bg  = torch.median(torch.concatenate(x))
else:
    mu_bg = torch.tensor(0, dtype=dtype, device=device)

########################################################

# EM for GMM parameters
print('Estimating GMM parameters')

print('  Resizing image, mask, and coordinates')
I_r, aff_r = my.torch_resize(Iim, aff, resolution, device, dtype=dtype)
M_r, _ = my.torch_resize(Mim, aff, resolution, device, dtype=dtype)
I_r[M_r<0.5] = mu_bg
Iim_shape = Iim.shape
del Iim, Mim, Sim, SimParc

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
# second bit: vox2vox transform to fake  space
T2 = my.getM(imageCOGvox[:-1,:], fakeCOGvox[:-1,:])
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

    normalizer = 1e-9 + torch.sum(W, dim=-1, keepdim=True)
    loglhood = torch.mean(torch.log(normalizer)).detach().cpu().numpy()
    W = W / normalizer

    # M step
    prior_loglhood = torch.zeros(1,dtype=dtype,device=device)
    for c in range(number_of_gmm_components.sum()):
        # crucially, we skip the background when we update the parameters (but we still add it to the cost)
        if c > 0:
            norm = torch.sum(W[:, :, :, c])
            means[c] = torch.sum(x * W[:, :, :, c]) / norm
            d = x - means[c]
            variances[c] = (torch.sum(d * d * W[:, :, :, c]) + prior_count * prior_variance) / ( norm + prior_count + 2)
        v = variances[c]
        prior_loglhood = prior_loglhood - ((1 + 0.5 * prior_count) * torch.log(v) + 0.5 * prior_count * prior_variance / v) / torch.numel(normalizer)
    loglhood = loglhood + prior_loglhood.detach().cpu().numpy()

    mixture_weights = torch.sum(W[:, :, :, 1:].reshape([np.prod(priors.shape[:-1]), number_of_gmm_components.sum() - 1]) + 1e-9, axis=0)
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
    -0.5 * torch.pow(I_r - means[c], 2.0) / variances[c]) + 1e-9

print('Computing normalizers (faster to do now with clustered priors)')
# We deform one class at the time; slower, but less memory
normalizers = torch.zeros(GAUSSIAN_LHOODS.shape[:-1], dtype=dtype, device=device)
# normalizers[h] = 1e-9
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
cmd = '  freeview -v ' + input_volume
if args.write_bias_corrected:
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






now2 = datetime.now()

current_time = now2.strftime("%H:%M:%S")
print("Current Time =", current_time)

runtime = now2 - now

print("Running Time =", runtime)


