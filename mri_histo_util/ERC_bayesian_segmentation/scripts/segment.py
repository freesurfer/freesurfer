#Imports
import os
import sys
basepath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(basepath)
from pathlib import Path
BASE_PATH = str(Path(__file__).resolve().parents[1])
sys.path.append(str(Path(__file__).resolve().parents[1]))
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
import ext.my_functions as my
from datetime import datetime
import ERC_bayesian_segmentation.relabeling as relab
import ERC_bayesian_segmentation.networks as nets
import glob
import scipy.sparse as sp
import ext.bias_field_correction_torch as bf
from ext.LBFGS import FullBatchLBFGS
import csv
import argparse
import ext.interpol as interpol
import math
import gc

########################################################

parser = argparse.ArgumentParser(description='Bayesian segmentation.')
parser.add_argument("--i", help="Image to segment.")
parser.add_argument("--i_seg", help="SynthSeg of image to segment (must include parcels). Will be computed if it does not exist")
parser.add_argument("--i_field", help="Registration to data/mni.nii.gz provided by mri_easyreg. Will be computed if if does not exist")
parser.add_argument("--atlas_dir", help="Atlas directory")
parser.add_argument("--hemi", help="hemisphere to segment (must be l or r)")
parser.add_argument("--gmm_mode", help="hemisphere to segment (must be l or r)", default="1mm")
parser.add_argument("--bf_mode", help="bias field basis function: dct, polynomial, or hybrid", default="dct")
parser.add_argument("--line_search", default="Armijo", help="type of line search; must be Armijo (default) or Wolfe")
parser.add_argument("--o", help="Output segmentation.")
parser.add_argument("--o_rgb", help="Output segmentation (soft, RGB).")
parser.add_argument("--o_atlas", help="Output deformed atlas (lumped classes)")
parser.add_argument("--o_vol", help="Output file with ROI volumes (csv format)")
parser.add_argument("--o_bf_corr", help="Bias field corrected volume, in .mgz or .nii(.gz) format")
parser.add_argument("--o_synthseg_vols", help="CSV file with SynthSeg volumes")
parser.add_argument("--cpu", action="store_true", help="Use CPU instead of GPU")
parser.add_argument("--threads", type=int, default=1, help="(optional) Number of CPU cores to be used. Default is 1. You can use -1 to use all available cores")
parser.add_argument("--skip_bf", action="store_true", help="Skip bias field correction")
parser.add_argument("--stiffness", type=float, default=4.5, help="(optional) Weight of regularized of atlas deformation")
args = parser.parse_args()

########################################################
if args.i is None:
    raise Exception('Input image is required')
if args.i_seg is None:
    raise Exception('SynthSeg file  of input image is required')
if args.i_field is None:
    raise Exception('MNI registration file of input image is required')
if args.atlas_dir is None:
    raise Exception('Atlas directory must be provided')
if args.hemi is None:
    raise Exception('Hemisphere must be provided')
if args.o is None:
    raise Exception('Output segmentation file must be provided')

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

# Input data
input_volume = args.i
input_seg = args.i_seg
mni_field = args.i_field
atlas_dir = args.atlas_dir
LUT_file = os.path.join(BASE_PATH, 'data', 'AllenAtlasLUT')
hemi = args.hemi
output_seg = args.o
output_seg_rgb = args.o_rgb
output_def_atlas = args.o_atlas
output_vol_file = args.o_vol
output_bf_corr = args.o_bf_corr
output_synthseg_vols = args.o_synthseg_vols
skip_bf = args.skip_bf
bf_mode = args.bf_mode
stiffness = args.stiffness

########################################################

# Constants
RESOLUTION_LEVELS = [3.0, 2.0, 1.0, 0.6666666666666666, 0.3333333333333333]
VOXEL_SKIP = [1, 1, 1, 1, 2]
NONLIN_CP_SPACING = [4, 4, 4, 6, 6]
dtype = torch.float32
SET_BG_TO_CSF = True # True = median of ventricles -> it seems much better than 0!

# The original numerical penalty and the analytical penalty match without scaling
# factors. The numerical penalty was divided by three and the analytical penalty
# is divided by two so the K need to be scaled with 3/2 to get the same regularization.
# Original K was 3 changed to 4.5. But it can be modified with --stiffness in the command line
K_PENALTY_GRAD = stiffness
RESOLUTION_ATLAS = 0.2
TOL = 1e-9
LR = 10.0
line_search = args.line_search
if line_search=='Armijo':
    STEPS = [1000, 800, 600, 400, 200] # we do EM every 25 of these
    # STEPS = [10, 10, 5] # Wolfe is slowerwe do EM every 25 of these
elif line_search == 'Wolfe':
    STEPS = [500, 400, 300, 200, 100] # Wolfe is slowerwe do EM every 25 of these
else:
    raise Exception('Line search must be Wolfe or Armijo')


if False:
    RESOLUTION_LEVELS = [3.0, 2.0, 1.0, 0.5]
    VOXEL_SKIP = [1, 1, 1, 2]
    NONLIN_CP_SPACING = [8, 8, 8, 8]
    STEPS = [200, 200, 200, 200]

if False:
    RESOLUTION_LEVELS = [3.0, 2.0, 1.0, 0.5, 0.2]
    VOXEL_SKIP = [1, 1, 1, 1, 1]
    NONLIN_CP_SPACING = [4, 4, 4, 6, 8]
    os.environ["CUDA_VISIBLE_DEVICES"] = ""
    STEPS = [1000, 800, 600, 400, 200]


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

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

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
    if (output_synthseg_vols is None) or (os.path.isfile(output_synthseg_vols)):
        synthseg_exists = True
    else:
        print('However, I did not find the CSV file with the volumes; I need to rerun SynthSeg, sorry!')

if synthseg_exists:
    print('Input segmentation exists; making sure it includes parcellation!')
    tmp, _ = my.MRIread(input_seg)
    if np.sum(tmp>1000)==0:
        raise Exception('Segmentation does not include parcellation! Please use different file or re-run SynthSeg with --parc')
else:
    print('Running SynthSeg')
    cmd = 'mri_synthseg --i ' + input_volume +  ' --o '  + input_seg + ' --threads ' + str(args.threads) + ' --cpu --parc --robust '
    if output_synthseg_vols is not None:
        cmd = cmd + ' --vol ' + output_synthseg_vols
    a = os.system(cmd + ' >/dev/null')
    if a > 0:
        print('Error in mri_synthseg; exitting...')
        sys.exit(1)

TMP_RESAMPLED = output_seg + '.tmp.resampled.mgz'
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
    Iim, _ = bf.correct_bias(Iim, Sim, maxit=100, penalty=0.1, order=4, device=device, dtype=dtype, basis=bf_mode)
    if output_bf_corr is not None:
        my.MRIwrite(Iim, aff, output_bf_corr)

print('Normalizing intensities')
Iim = Iim * 110 / np.median(Iim[(Sim==2) | (Sim==41)])

########################################################

print('Subdividing brainstem into left and right halves')
BRAINSTEM = 16
BRAINSTEM_L = 161
BRAINSTEM_R = 162

if os.path.isfile(mni_field):
    print('   Registration file found; no need to run EasyReg!')
else:
    print('   Running EasyReg')
    mni_image = os.path.join(BASE_PATH, 'data', 'mni.nii.gz')
    mni_seg = os.path.join(BASE_PATH, 'data', 'mni.synthseg.nii.gz')
    a = os.system('mri_easyreg --flo ' + mni_image + ' --ref ' + input_volume +
                  ' --flo_seg ' + mni_seg + ' --ref_seg ' + input_seg + ' --threads ' + str(args.threads) + ' --fwd_field ' + mni_field + '  >/dev/null')
    if a>0:
        print('Error in mri_easyreg; exitting...')
        sys.exit(1)
FIELD, _ = my.MRIread(mni_field)
LEFT = (FIELD[:,:,:,0]<0)
Sim[(Sim==BRAINSTEM) & LEFT] = BRAINSTEM_L
Sim[(Sim==BRAINSTEM) & (LEFT==0)] = BRAINSTEM_R


########################################################

print('Creating and applying mask')
if hemi=='l':
    M = ( (Sim < 30) | (Sim==161) | ( (Sim > 1000) & (Sim < 2000) ) )
elif hemi=='r':
    M = ( ( (Sim > 40) & (Sim<100) )  | (Sim==162) | (Sim>2000) )
if True:
    M[Sim==4] = 0
    M[Sim==5] = 0
    M[Sim==43] = 0
    M[Sim==44] = 0
M[Sim==14] = 0
M[Sim==15] = 0
M[Sim==24] = 0
M[Sim==0] = 0
# I now do this with the resampled mask at each level, to avoid blurring mask edges
# Iim[M==0] = 0
Mim, cropping = my.cropLabelVol(M, margin=5)
Sim = my.applyCropping(Sim, cropping)
Iim = my.applyCropping(Iim, cropping)
aff[:3, -1] = aff[:3, -1] + aff[:-1, :-1] @ cropping[:3]

########################################################

# Read atlas and merge classes
print('Reading in atlas')

# Get the label groupings and atlas labels
# from the config files
aseg_label_list = np.unique(Sim[Mim>0])
tissue_index, grouping_labels, label_list, number_of_gmm_components = relab.get_tissue_settings(
            os.path.join(BASE_PATH, 'data', 'atlas_names_and_labels.yaml'),
            os.path.join(BASE_PATH, 'data', 'combined_atlas_labels_' + args.gmm_mode + '.yaml'),
            os.path.join(BASE_PATH, 'data', 'combined_aseg_labels_' + args.gmm_mode + '.yaml'),
            os.path.join(BASE_PATH, 'data', 'gmm_components_' + args.gmm_mode + '.yaml'),
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

aff_A = np.diag([.2, .2, .2, 1])
if hemi=='r':
    aff_A = np.diag([-.2, .2, .2, 1])

########################################################
print('Computing initial values for means and variances')
mus_ini = []
vars_ini = []
mixture_weights = []

for t in range(len(number_of_gmm_components)-1):
    labs = grouping_labels[t+1]
    x = []
    for l in labs:
        x.append(Iim[Sim==l])

    if len(x) > 0:
        x = np.concatenate(x)
        mu = np.median(x)
        std = 1.4826 * np.median(np.abs(x - mu))
        var = std ** 2
        if number_of_gmm_components[t+1]==1:
            mus_ini.append(mu[np.newaxis])
            vars_ini.append(var[np.newaxis])
            mixture_weights.append(np.ones(1))
        else:
            # Estimate GMM with shared variance (avoids a component with tiny variance)
            nc = number_of_gmm_components[t+1]
            nx = len(x)
            gmm_mus = np.linspace(mu - 0.5 * std, mu + 0.5 * std, nc)
            gmm_var= var * np.ones(1)
            gmm_ws = (1 / float(nc)) * np.ones(nc)
            W = np.zeros([nx, nc])
            for its in range(200):
                # E step
                for c in range(nc):
                    W[:, c] = gmm_ws[c] / np.sqrt(2.0 * np.pi * np.sqrt(gmm_var)) * np.exp(-0.5 * (x - gmm_mus[c])**2 / gmm_var)
                normalizer = np.sum(W + 1e-9, axis=1)
                # print(-np.mean(np.log(normalizer)))
                W /= normalizer[:, np.newaxis]
                # M step
                denominators = np.sum(W, axis=0)
                gmm_ws = denominators / np.sum(denominators)
                gmm_var = 0
                for c in range(nc):
                    gmm_mus[c] = np.sum(W[:, c] * x) / denominators[c]
                    aux = x - gmm_mus[c]
                    gmm_var += np.sum(W[:, c] * aux * aux)
                gmm_var /= np.sum(denominators)

            mus_ini.append(gmm_mus)
            vars_ini.append(gmm_var * np.ones(nc))
            mixture_weights.append(gmm_ws)

mus_ini = np.concatenate(mus_ini)
vars_ini = np.concatenate(vars_ini)
mixture_weights = np.concatenate(mixture_weights)

if SET_BG_TO_CSF:
    x = []
    for l in [4, 5, 43, 44]: # , 24]:
        x.append(Iim[Sim==l])
    mu_bg  = np.median(np.concatenate(x))
else:
    mu_bg = 0

########################################################

print('Coarse alignment by matching centers of gravity')
idx = np.where(A[:,:,:,0]<0.5)
cog_atl_vox = np.array([[np.mean(idx[0])], [np.mean(idx[1])], [np.mean(idx[2])]])
cog_atl_ras = np.squeeze(my.vox2ras(cog_atl_vox, aff_A))

idx = np.where(Mim>0)
cog_mri_vox = np.array([[np.mean(idx[0])], [np.mean(idx[1])], [np.mean(idx[2])]])
cog_mri_ras = np.squeeze(my.vox2ras(cog_mri_vox, aff))

aff_A[:-1,-1] = aff_A[:-1,-1] - np.squeeze(cog_atl_ras)
aff[:-1,-1] = aff[:-1,-1] - np.squeeze(cog_mri_ras)


########################################################

# OK so we work one resolution at the time
n_levels = len(RESOLUTION_LEVELS)

# Initialize
ts = None
thetas = None
shears = None
scalings = None
mus = mus_ini
vars = vars_ini
var_bg = np.min(vars_ini)
FIELD = None
model = None

for n in range(n_levels):

    print('*******************')
    print('*   Level ' + str(n+1) + ' / ' + str(n_levels) + '   *')
    print('*******************')

    print('Resizing image')
    I_r, aff_r = my.torch_resize(Iim, aff, RESOLUTION_LEVELS[n], device, dtype=dtype)
    print('Resizing mask')
    M_r, _ = my.torch_resize(Mim, aff, RESOLUTION_LEVELS[n], device, dtype=dtype)
    I_r[M_r<0.5] = mu_bg
    print('Resizing atlas')
    slow = (device.type == 'cpu' or torch.cuda.get_device_properties(device).total_memory < 30*1024**3)
    A_r, aff_A_r = my.torch_resize(A, aff_A, RESOLUTION_LEVELS[n], device, dtype=dtype, slow=slow)

    # Go over resolutions / modes
    for flexibility_mode in range(3):

        if flexibility_mode == 0:
            allow_similarity = True
            allow_shearing = False
            allow_nonlinear = False
            allow_gaussian_update = False
            print('First pass: similarity mode')
        elif flexibility_mode == 1:
            allow_similarity = True
            allow_shearing = True
            allow_nonlinear = False
            allow_gaussian_update = False
            print('Second pass: affine')
        else:
            allow_similarity = False
            allow_shearing = False
            # Setting the nonlinear and gaussian updates to same condition
            allow_nonlinear = (RESOLUTION_LEVELS[n]<=2.25)
            allow_gaussian_update = (RESOLUTION_LEVELS[n]<=1.25)
            print('Third pass: affine + nonlinear')

        if allow_gaussian_update:
            outter_its = int(STEPS[n]/25)
            inner_its = 25
        else:
            outter_its = 1
            inner_its = STEPS[n]
        if (allow_similarity==False) and (allow_shearing==False) and (allow_nonlinear==False):
            outter_its = 0 # if nothing to do, just skip

        # Outter loop over deformation / EM
        for outter_it in range(outter_its):

            # Compute cost of Gaussian parameters
            prior_count = 100 / (RESOLUTION_LEVELS[n] ** 3) / (VOXEL_SKIP[n] ** 3)
            prior_variance = var_bg
            denominator = np.prod(I_r[::VOXEL_SKIP[n], ::VOXEL_SKIP[n], ::VOXEL_SKIP[n]].shape)
            prior_loglhood = - ((1 + 0.5 * prior_count) * np.log(var_bg) + 0.5 * prior_count * prior_variance / var_bg) / denominator
            for c in range(len(vars)):
                prior_loglhood -= ((1 + 0.5 * prior_count) * np.log(vars[c]) + 0.5 * prior_count * prior_variance / vars[c]) / denominator

            print('   Deforming atlas: outter iteration ' + str(outter_it+1) + ' of ' + str(outter_its))

            model = nets.AtlasDeformer(I_r, aff_r, A_r, aff_A_r, ts_ini=ts,
                                        thetas_ini=thetas, scalings_ini=scalings,
                                        shears_ini=shears, mus=mus, vars=vars,
                                        mu_bg=mu_bg, var_bg=var_bg, weights=mixture_weights,
                                        gmm_components=number_of_gmm_components[1:],
                                        k_penalty_grad=K_PENALTY_GRAD,
                                        FIELD_ini=FIELD,
                                        cp_spacing=NONLIN_CP_SPACING[n],
                                        allow_similarity=allow_similarity,
                                        allow_shearing=allow_shearing,
                                        allow_nonlinear=allow_nonlinear,device=device,
                                        skip=VOXEL_SKIP[n], dtype=dtype)

            optimizer = FullBatchLBFGS(model.parameters(), lr=LR, line_search=line_search)

            loss = model()[0]
            loss.backward()

            loss_old = 1e10
            for inner_it in range(inner_its):

                # define closure for line search
                def closure():
                    optimizer.zero_grad()
                    loss_fn = model()[0]
                    return loss_fn


                # perform line search step
                options = {'closure': closure, 'current_loss': loss, 'eta': 2, 'max_ls': 100,
                            'interpolate': True, 'inplace': False}
                loss = (optimizer.step(options=options))[0]
                # compute gradient at new iterate (only for Armijo)
                if line_search=='Armijo':
                    loss.backward()

                # print step info
                loss_np = loss.detach().cpu().numpy() - prior_loglhood
                print('         Step %d of atlas deformation, loss = %.6f' % (inner_it + 1, loss_np), flush=True)

                if ((loss_old - loss_np) < TOL):
                    print('         Decrease in loss below tolerance limit')
                    break
                else:
                    loss_old = loss_np

            # Retrieve model parameters
            ts = model.ts.detach().cpu().numpy()
            thetas = model.thetas.detach().cpu().numpy()
            shears = model.shears.detach().cpu().numpy()
            scalings = model.scalings.detach().cpu().numpy()
            FIELD = model.FIELD.detach().cpu().numpy()
            loss = loss_np
            if allow_gaussian_update:
                # priors = model()[2].detach().cpu().numpy()
                _, _, priors, deformation_penalty = model()

            del model
            del optimizer
            gc.collect()
            torch.cuda.empty_cache()

            # Now we do EM if needed
            if allow_gaussian_update:
                with torch.no_grad():
                    print('   Updating Gaussian parameters with EM: outter iteration ' + str(outter_it+1) + ' of ' + str(outter_its))
                    means = torch.tensor([mu_bg, *mus], device=device, dtype=dtype)
                    variances = torch.tensor([var_bg, *vars], device=device, dtype=dtype)
                    weights = torch.tensor(np.array([1.0] + mixture_weights.tolist()), device=device, dtype=dtype)
                    priors = torch.tensor(priors, device=device, dtype=dtype)
                    W = torch.zeros([*priors.shape[:-1], number_of_gmm_components.sum()]).to(device)
                    x = torch.tensor(I_r[::VOXEL_SKIP[n], ::VOXEL_SKIP[n], ::VOXEL_SKIP[n]], device=device, dtype=dtype)

                    # We now put a Scaled inverse chi-squared prior on the variances to prevent them from going to zero
                    prior_count = 100 / (RESOLUTION_LEVELS[n] ** 3) / (VOXEL_SKIP[n] ** 3)
                    prior_variance =var_bg
                    loglhood_old = -10000
                    for em_it in range(10):
                        # E step
                        for c in range(n_tissues):
                            prior = priors[:,:,:,c]
                            num_components = number_of_gmm_components[c]
                            for g in range(num_components):
                                gaussian_number = sum(number_of_gmm_components[:c]) + g
                                d = x - means[gaussian_number]
                                W[:,:,:,gaussian_number] = weights[gaussian_number] * prior * torch.exp(-d * d / (2 * variances[gaussian_number])) / torch.sqrt(2.0*math.pi*variances[gaussian_number])

                        normalizer = 1e-9 + torch.sum(W, dim=-1, keepdim=True)
                        loglhood = torch.mean(torch.log(normalizer)).detach().cpu().numpy()
                        W = W / normalizer

                        # M step
                        prior_loglhood = np.zeros_like(loglhood)
                        for c in range(number_of_gmm_components.sum()):
                            # crucially, we skip the background when we update the parameters (but we still add it to the cost)
                            if c>0:
                                norm = torch.sum(W[:,:,:,c])
                                means[c] = torch.sum(x * W[:,:,:,c]) / norm
                                d = x - means[c]
                                variances[c] = (torch.sum(d * d * W[:, :, :, c])  + prior_count * prior_variance ) / (norm + prior_count + 2)
                            v = variances[c].detach().cpu().numpy()
                            prior_loglhood = prior_loglhood - (  (1 + 0.5 * prior_count) *  np.log(v) +   0.5 * prior_count * prior_variance / v   )/ torch.numel(normalizer)
                        loglhood = loglhood + prior_loglhood

                        mixture_weights = torch.sum(W[:, :, :, 1:].reshape([np.prod(priors.shape[:-1]), number_of_gmm_components.sum()-1]) + 1e-9, axis=0)
                        for c in range(n_tissues-1):
                            # mixture weights are normalized (those belonging to one mixture sum to one)
                            num_components = number_of_gmm_components[c+1]
                            gaussian_numbers = torch.tensor(np.sum(number_of_gmm_components[1:c+1]) + \
                                                            np.array(range(num_components)), device=device, dtype=dtype).int()

                            mixture_weights[gaussian_numbers] /= torch.sum(mixture_weights[gaussian_numbers])

                        weights[1:] = mixture_weights

                        if (torch.sum(torch.isnan(means))>0) or (torch.sum(torch.isnan(variances))>0):
                            print('nan in Gaussian parameters...')
                            import pdb; pdb.set_trace()

                        print('         Step %d of EM, -loglhood = %.6f' % (em_it + 1, -loglhood + deformation_penalty ), flush=True)
                        if (loglhood-loglhood_old)<TOL:
                            print('         Decrease in loss below tolerance limit')
                            break
                        else:
                            loglhood_old = loglhood
                    mus = means.detach().cpu().numpy()[1:]
                    vars = variances.detach().cpu().numpy()[1:]
                    del W
                    del x
                    del priors
                    gc.collect()
                torch.cuda.empty_cache()


########################################################

# Compute final outputs at finest resolution with the latest model
with torch.no_grad():

    # TODO: without this line, I get weird runtime errors...
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    # At this point we need just the likelihood and the normalizer. The lumped priors are
    # used just for debugging/checking so I'll keep them for now, but only on the CPU.
    # In any case, we don't need to run the whole model again, so let's abuse the constructor
    # to get the necessary stuff prepared and then only get the grids, normalizer, and usampled priors
    model = nets.AtlasDeformer(I_r, aff_r, A_r, aff_A_r, ts_ini=ts, thetas_ini=thetas, scalings_ini=scalings,
                               shears_ini=shears, mus=mus, vars=vars, mu_bg=mu_bg, var_bg=var_bg,
                               weights=mixture_weights, gmm_components=number_of_gmm_components[1:],
                               k_penalty_grad=K_PENALTY_GRAD,
                               FIELD_ini=FIELD, cp_spacing=NONLIN_CP_SPACING[n], allow_similarity=False,
                               allow_shearing=False, allow_nonlinear=False, device=device, dtype=dtype)
    print('Deforming grouped atlas and computing normalizer')
    GAUSSIAN_LHOODS = model.GAUSSIAN_LHOODS
    grids, _ = model.get_atlas_sampling_grid()
    normalizer, priors = model.get_priors(grids, return_llhood=False)
    del model
    gc.collect()
    torch.cuda.empty_cache()

    # my.viewVolume([I_r, priors], aff_r)

    print('Deforming one label at the time')
    seg = torch.zeros(normalizer.shape, dtype=torch.int, device=device)
    seg_rgb = torch.zeros([*normalizer.shape, 3], dtype=dtype, device=device)
    max_p = torch.zeros(normalizer.shape, dtype=dtype, device=device)
    vols = torch.zeros(n_labels, device=device, dtype=dtype)
    names, colors = my.read_LUT(LUT_file)
    voxel_vol = np.abs(np.linalg.det(aff_r))

    # TODO: choose good number of workers/prefetch factor
    for n, (prior_indices, prior_values) in enumerate(label_loader):
        print('Reading in label ' + str(n + 1) + ' of ' + str(n_labels), end='\r', flush=True)

        if prior_indices.numel() == 0:
            continue
        prior_indices = torch.as_tensor(prior_indices, device=device, dtype=torch.long).squeeze()
        prior_values = torch.as_tensor(prior_values, device=device, dtype=dtype).squeeze()

        aff_pr = np.copy(aff_A)
        if n == 0:
            # background
            prior = torch.sparse_coo_tensor(prior_indices[None], prior_values,
                                            [torch.Size(atlas_size).numel()]).to_dense()
            del prior_indices, prior_values
            prior = prior.reshape(torch.Size(atlas_size))
            lr_crop = (slice(None),) * 3
        else:
            # find bounding box of label in atlas space
            prior_indices = my.ind2sub(prior_indices, atlas_size)
            min_x, max_x = prior_indices[0].min().item(), prior_indices[0].max().item() + 1
            min_y, max_y = prior_indices[1].min().item(), prior_indices[1].max().item() + 1
            min_z, max_z = prior_indices[2].min().item(), prior_indices[2].max().item() + 1
            crop_atlas_size = [max_x - min_x, max_y - min_y, max_z - min_z]
            prior_indices[0] -= min_x
            prior_indices[1] -= min_y
            prior_indices[2] -= min_z
            prior = torch.sparse_coo_tensor(prior_indices, prior_values, crop_atlas_size).to_dense()
            del prior_indices, prior_values
            aff_pr[:3, -1] += aff_pr[:3, :3] @ np.asarray([min_x, min_y, min_z])

            # find bounding box of label in MRI space
            hr2lr = np.linalg.inv(aff_A_r) @ aff_A
            min_x, min_y, min_z = (hr2lr[:3, :3] @ np.asarray([min_x-1, min_y-1, min_z-1] + hr2lr[:3, -1])).tolist()
            max_x, max_y, max_z = (hr2lr[:3, :3] @ np.asarray([max_x, max_y, max_z] + hr2lr[:3, -1])).tolist()
            mask =  (grids[0, ..., 0] >= min_x)
            mask &= (grids[0, ..., 0] <= max_x)
            mask &= (grids[0, ..., 1] >= min_y)
            mask &= (grids[0, ..., 1] <= max_y)
            mask &= (grids[0, ..., 2] >= min_z)
            mask &= (grids[0, ..., 2] <= max_z)
            if ~mask.any():
                continue
            nx, ny, nz = mask.shape
            tmp = mask.reshape([nx, -1]).any(-1).nonzero()
            lr_min_x, lr_max_x = tmp.min().item(), tmp.max().item() + 1
            tmp = mask.movedim(0, -1).reshape([ny, -1]).any(-1).nonzero()
            lr_min_y, lr_max_y = tmp.min().item(), tmp.max().item() + 1
            tmp = mask.reshape([-1, nz]).any(0).nonzero()
            lr_min_z, lr_max_z = tmp.min().item(), tmp.max().item() + 1
            del tmp, mask
            lr_crop = (slice(lr_min_x, lr_max_x), slice(lr_min_y, lr_max_y), slice(lr_min_z, lr_max_z))

        if RESOLUTION_LEVELS[-1] != RESOLUTION_ATLAS:
            prior, aff_pr = my.torch_resize(prior, aff_pr, RESOLUTION_LEVELS[-1], device, dtype=dtype)

        # shift/scale sampling grid appropriately
        aff_shift = torch.as_tensor(np.linalg.inv(aff_pr) @ aff_A_r, device=device, dtype=dtype)
        grids_shifted = grids[(slice(None), *lr_crop, slice(None))]
        grids_shifted = aff_shift[:3, :3].matmul(grids_shifted.unsqueeze(-1)).squeeze(-1)
        grids_shifted = grids_shifted.add_(aff_shift[:3, -1])

        prior = interpol.grid_pull(prior[None, None, :, :, :], grids_shifted, interpolation=1)
        del grids_shifted

        num_components = number_of_gmm_components[tissue_index[n]]
        gaussian_numbers = torch.tensor(np.sum(number_of_gmm_components[:tissue_index[n]]) + \
                                        np.array(range(num_components)), device=device, dtype=dtype).int()
        lhood = torch.sum(GAUSSIAN_LHOODS[:, :, :, gaussian_numbers] * weights[None, None, None, gaussian_numbers], 3)
        post = torch.squeeze(prior)
        post *= lhood[lr_crop]
        post /= normalizer[lr_crop]
        del prior

        vols[n] = torch.sum(post) * voxel_vol
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

# Remeber to undo the centering of the center of mass!
aff[:-1,-1] = aff[:-1,-1] + np.squeeze(cog_mri_ras)
affine = aff_r.copy()
affine[:-1, -1] = affine[:-1, -1] + np.squeeze(cog_mri_ras)

my.MRIwrite(seg.detach().cpu().numpy(), affine, output_seg, dtype='int')
if output_seg_rgb is not None:
    my.MRIwrite(seg_rgb.detach().cpu().numpy(), affine, output_seg_rgb)
if output_def_atlas is not None:
    my.MRIwrite(priors, affine, output_def_atlas)
if output_vol_file is not None:
    vols = vols.detach().cpu().numpy()
    with open(output_vol_file, 'w') as csvfile:
        writer = csv.writer(csvfile)
        aux = label_list[1:]
        row = []
        for l in aux:
            row.append(names[int(l)])
        writer.writerow(row)
        row = []
        for j in range(1, 334):
            row.append(str(vols[j]))
        writer.writerow(row)

cmd = 'freeview -v ' + input_volume
if output_bf_corr is not None:
    cmd = cmd + ' -v ' + output_bf_corr
cmd = cmd + ' -v ' + output_seg + ':colormap=lut:lut=' + LUT_file
if output_seg_rgb is not None:
    cmd = cmd + ' -v ' + output_seg_rgb + ':rgb=true'
if output_def_atlas is not None:
    cmd = cmd + ' -v ' + output_def_atlas + ':colormap=heat'
print(cmd)
if output_vol_file is not None:
    print('oocalc %s' % (output_vol_file))

print('All done!')

now2 = datetime.now()

current_time = now2.strftime("%H:%M:%S")
print("Current Time =", current_time)

runtime = now2 - now

print("Running Time =", runtime)




