from argparse import ArgumentParser
import os
import pdb
from collections import OrderedDict

import numpy as np
import nibabel as nib
from scipy.special import softmax
from scipy.ndimage import distance_transform_edt, gaussian_filter, binary_dilation, generate_binary_structure
import torch

# Need it if providing the posteriors.
ASEG_LABELS= {
    'Background': 0,
    'Right-Hippocampus': 53,
    'Left-Hippocampus': 17,
    'Right-Lateral-Ventricle': 43,
    'Left-Lateral-Ventricle': 4,
    'Right-Thalamus': 49,
    'Left-Thalamus': 10,
    'Right-Amygdala': 54,
    'Left-Amygdala': 18,
    'Right-Putamen': 51,
    'Left-Putamen': 12,
    'Right-Pallidum': 52,
    'Left-Pallidum': 13,
    'Right-Cerebrum-WM': 41,
    'Left-Cerebrum-WM': 2,
    'Right-Cerebellar-WM': 46,
    'Left-Cerebellar-WM': 7,
    'Right-Cerebrum-GM': 42,
    'Left-Cerebrum-GM': 3,
    'Right-Cerebellar-GM': 47,
    'Left-Cerebellar-GM': 8,
    'Right-Caudate': 50,
    'Left-Caudate': 11,
    'Brainstem': 16,
    '4th-Ventricle': 15,
    '3rd-Ventricle': 14,
    'Right-Accumbens': 58,
    'Left-Accumbens': 26,
    'Right-VentralDC': 60,
    'Left-VentralDC': 28,
    'Right-Inf-Lat-Ventricle': 44,
    'Left-Inf-Lat-Ventricle': 5,
}

# Need it for clustering regions under the same gaussian.
CLUSTER_DICT = {
    'Gray': [53, 17, 51, 12, 54, 18, 50, 11, 58, 26, 42, 3],
    'CSF': [4, 5, 43, 44, 15, 14, 24],
    'Thalaumus': [49, 10],
    'Pallidum': [52, 13],
    'VentralDC': [28, 60],
    'Brainstem': [16],
    'WM': [41, 2],
    'cllGM': [47, 8],
    'cllWM': [46, 7]
}

# make polynomial basis functions
def get_basis_functions(shape, order=3, device='cpu', dtype=torch.float32):

    G = torch.meshgrid(torch.arange(shape[0]), torch.arange(shape[1]), torch.arange(shape[2]))
    Gnorm = []
    for i in range(3):
        aux = G[i].type(dtype).to(device)
        Gnorm.append(2 * ((aux / (shape[i] - 1)) - 0.5))

    B = []

    for x in range(order + 1):
        for y in range(order + 1):
            for z in range(order + 1):
                if ((x + y + z) <= order) and ((x + y + z) > 0):
                    b = torch.ones(shape, device=device, dtype=dtype)
                    for i in range(x):
                        b = b * Gnorm[0]
                    for i in range(y):
                        b = b * Gnorm[1]
                    for i in range(z):
                        b = b * Gnorm[2]
                    B.append(b)
    return B

# make polynomial basis functions
def get_basis_functions_dct(shape, order=3, device='cpu', dtype=torch.float32):

    one_d_basis_x = []
    one_d_basis_y = []
    one_d_basis_z = []
    for i in range(order):
        one_d_basis_x.append(torch.tensor(np.cos((2.0 * np.arange(shape[0]) + 1) * np.pi * (i + 1) / (2.0 * shape[0])), device=device, dtype=dtype))
        one_d_basis_y.append(torch.tensor(np.cos((2.0 * np.arange(shape[1]) + 1) * np.pi * (i + 1) / (2.0 * shape[1])), device=device, dtype=dtype))
        one_d_basis_z.append(torch.tensor(np.cos((2.0 * np.arange(shape[2]) + 1) * np.pi * (i + 1) / (2.0 * shape[2])), device=device, dtype=dtype))

    B = []
    for x in range(order + 1):
        for y in range(order + 1):
            for z in range(order + 1):
                if ((x + y + z) <= order) and ((x + y + z) > 0):
                    b = torch.ones(shape, device=device, dtype=dtype)
                    for i in range(x):
                        b *= one_d_basis_x[i][:, None, None]
                    for i in range(y):
                        b *= one_d_basis_y[i][None, :, None]
                    for i in range(z):
                        b *= one_d_basis_z[i][None, None, :]
                    B.append(b)
    return B

# Main function to correct bias field
def correct_bias(mri, seg, maxit=100, penalty=0.1, order=4, basis='hybrid', device='cpu', dtype=torch.float32):

    with torch.no_grad():
        # get image and masks as tensors, all masked by segmentation
        mri[mri<0] = 0 # avoid logarithm of negative numbers...
        mask = seg > 0
        I = torch.tensor(np.squeeze(mri)[mask > 0], device=device, dtype=dtype)
        nvox = I.shape[0]
        nclass = len(CLUSTER_DICT)
        post = torch.zeros([nvox, nclass], device=device, dtype=dtype)
        seg[seg >= 2000] = 42
        seg[seg > 1000] = 3

        print('Gaussian filtering for bias field correction')
        sigma = 1.0
        sl = np.ceil(sigma * 2.5).astype(int)
        v = np.arange(-sl, sl + 1)
        gauss = np.exp((-(v / sigma) ** 2 / 2))
        kernel = gauss / np.sum(gauss)
        kernel = torch.tensor(kernel, device=device, dtype=dtype)
        kernel = kernel[None, None, None, None, :]

        for it_lab, (lab_str, lab_list) in enumerate(CLUSTER_DICT.items()):
            M = torch.zeros(mri.shape, device=device, dtype=dtype)
            for lab in lab_list:
                M[seg==lab] = 1.0
            M = M[None, None, :,:,:]
            for d in range(3):
                M = torch.conv3d(M, kernel, bias=None, stride=1, padding=[0, 0, int((kernel.shape[-1] - 1) / 2)])
                M = M.permute([0, 1, 4, 2, 3])
            M = torch.squeeze(M)
            post[:, it_lab] = M[mask]

        # Get basis functions, and mask by segmentation as well
        if basis=='hybrid':
            print('  Using hybrid (DCT+polynomial) basis functions')
            BFs = get_basis_functions_dct(mask.shape, order, device, dtype=dtype) + get_basis_functions(mask.shape, order, device, dtype=dtype)
        elif basis=='dct':
            print('  Using DCT basis functions')
            BFs = get_basis_functions_dct(mask.shape, order, device, dtype=dtype)
        elif basis == 'polynomial':
            print('  Using polynomial basis functions')
            BFs = get_basis_functions(mask.shape, order, device, dtype=dtype)
        else:
            raise Exception('basis must be dct, polynomial, or hybrid')
        nbf = len(BFs)
        A = torch.zeros([nvox, nbf], device=device, dtype=dtype)
        for i in range(nbf):
            A[:, i] = BFs[i][mask]

        # Log transform with scaling
        factor = 1000.0 / torch.max(I)
        y = torch.log(1 + I * factor)

        # Main loop
        print('Bias field correction')
        C = torch.zeros(nbf, device=device, dtype=dtype)
        wij = torch.zeros([nvox, nclass], device=device, dtype=dtype)
        R = torch.zeros(nvox, device=device, dtype=dtype)
        REG = penalty * torch.eye(nbf, device=device, dtype=dtype)
        mus = torch.zeros(nclass, device=device, dtype=dtype)
        vars = torch.zeros(nclass, device=device, dtype=dtype)
        normalizer = torch.sum(post, axis=0)
        ready = False
        it = 0
        while ready==False:
            it = it + 1
            ycorr = y - torch.sum(A * C, dim=1)
            cost = torch.zeros(1, device=device, dtype=dtype)
            for j in range(nclass):
                mus[j]  = torch.sum(ycorr * post[:, j]) / normalizer[j]
                aux = ycorr - mus[j]
                vars[j] = torch.sum( (aux * aux) * post[:, j]) / normalizer[j]
                wij[:, j] = post[:, j] / vars[j]
                cost += torch.mean(post[:,j] * (0.5 * torch.log(2*torch.pi*vars[j]) +0.5 * (ycorr - mus[j])**2 ))
            wi = torch.sum(wij, axis=1)

            R[:] = y[:]
            for j in range(nclass):
                R -=  ((wij[:,j] / wi) * mus[j])

            Cold = torch.clone(C)
            C = torch.inverse(A.T @ (wi[..., None] * A) + REG) @ (A.T @ (wi * R))
            diff = torch.sum((C-Cold)**2)
            print('  Iteration ' + str(it) + ': difference is ' + str(diff.item()))
            if diff<1e-6:
                print('  Converged')
                ready = True
            if it==maxit:
                print('  Tired convergence')
                ready = True

        Icorr = torch.zeros(mask.shape, device=device, dtype=dtype)
        Icorr[mask] = (torch.exp(ycorr) - 1) / factor

        Icorr = Icorr.detach().cpu().numpy()
        cost = cost.detach().cpu().numpy()

    torch.cuda.empty_cache()

    return Icorr, cost
