import torch
import numpy as np
from torch.nn.functional import conv3d
from  ext.unet3d.model import UNet2D

def photo_imputation(input, aff, y_shifts, sz, threshold_fg, unsharp_sigma, unsharp_amount, checkpoint, device):
    with torch.no_grad():

        # Prepare CNN
        model = UNet2D(4, 1, final_sigmoid=False, f_maps=128, layer_order='gcl', num_groups=8, num_levels=5,
                       is_segmentation=False).to(device)
        if device.type == 'cpu':
            cp = torch.load(checkpoint, map_location=torch.device('cpu'))
        else:
            cp = torch.load(checkpoint)
        model.load_state_dict(cp['model_state_dict'])

        # Prepare normalized input tensor
        I = input.clone()
        affine = aff.copy()
        volres = np.sqrt(np.sum(aff**2, axis=0))[:-1]
        if np.any(volres[:2]!=1.0):
            resized = []
            from cv2 import resize, INTER_AREA
            for k in range(I.shape[2]):
                Isl = resize(I[:,:,k,:].detach().cpu().numpy(), None, fx=volres[0], fy=volres[1], interpolation=INTER_AREA)
                resized.append(Isl)
            I = torch.tensor(np.stack(resized, axis=-2), device=device, dtype=I.dtype)
            affine[:-1, 0] /= volres[0]
            affine[:-1, 1] /= volres[1]
            affine[:-1, -1] += (affine[:-1, :-1] @ np.array([0.5-0.5*volres[0], 0.5-0.5*volres[1], 0]))

        Igray = I.mean(dim=3)
        M = (Igray > threshold_fg)

        # normalize slices by their medians
        medians_slices = torch.zeros(Igray.shape[2], dtype=torch.float, device=device)
        for j in range(Igray.shape[2]):
            im = Igray[:, :, j]
            ma = M[:, :, j]
            if ma.sum() > 0:
                medians_slices[j] = torch.median(im[ma])
        medians_slices = medians_slices / torch.mean(medians_slices[medians_slices > 0])
        medians_slices[medians_slices == 0] = 1
        I /= medians_slices[None, None, :, None]

        # normalize input by median of whole thing (per channel)
        medians_whole = torch.zeros(3, device=device, dtype=I.dtype)
        for c in range(3):
            auxI = I[:, :, :, c]
            medians_whole[c] = torch.median(auxI[M])
            I[:, :, :, c] /= medians_whole[c]

        # Allocate input tensor with 3 batches (rgb), 4 channels (sliceA, sliceB, distA, distB), and 32x size
        shape2d = np.array(I.shape[:2])
        W = (np.ceil(shape2d / 32.0) * 32).astype('int')
        idx = np.floor((W - shape2d) / 2).astype('int')
        S = torch.zeros([3, 4, *W], dtype=torch.float32, device=device)  # will store inputs

        # Prepare list of new coordinates
        affnew = affine.copy()
        affnew[1, 2] /= np.abs(affnew[1, 2])
        affnew[1, 3] = affnew[1, 3] + 0.5 * (np.abs(affine[1, 2]) - 1)
        nk = np.ceil(np.abs(affine[1, 2]) * I.shape[2]).astype(np.int32)
        yp = affine[1, 3] + affine[1, 2] * np.arange(I.shape[2])
        if y_shifts is not None:
            yp -= (y_shifts * np.exp(sz / 20))  # crucial to scale shifts

        # Allocate outputs
        LINEAR = torch.zeros([*I.shape[:2], nk, 3], dtype=torch.float32, device=device)
        PRED = torch.zeros_like(LINEAR)

        # Loop over slices to interpolate!
        for k in range(nk):
            print('Working on slice ' + str(k + 1) + ' of ' + str(nk), end='\r')
            y = affnew[1, 2] * k + affnew[1, 3]
            idx1 = np.where(yp >= y)[0]
            if (len(idx1) > 0) and (len(idx1) < len(yp)):  # otherwise nothing to do
                idx1 = idx1.max()
                idx2 = min(idx1 + 1, len(yp) - 1)
                d1 = (yp[idx1] - y).clip(1e-6)
                d2 = (y - yp[idx2]).clip(1e-6)
                w1 = d2 / (d1 + d2)
                w2 = 1.0 - w1
                linear_interp = w1 * I[:, :, idx1, :] + w2 * I[:, :, idx2, :]
                LINEAR[:, :, k, :] = linear_interp
                for c in range(3):  # note that we trained with slices from back to front therefore the reordering
                    S[c, 1, idx[0]:idx[0] + I.shape[0], idx[1]:idx[1] + I.shape[1]] = I[:, :, idx1, c]
                    S[c, 0, idx[0]:idx[0] + I.shape[0], idx[1]:idx[1] + I.shape[1]] = I[:, :, idx2, c]
                    S[c, 3, :, :] = 0.1 * d1
                    S[c, 2, :, :] = 0.1 * d2
                pred = model(S.flip(dims=(2, 3)).permute(0, 1, 3, 2)).flip(dims=(2, 3)).permute(0, 1, 3, 2)  # again, because we trained with identity vox2ras
                PRED[:, :, k, :] = pred[:, 0, :, :].permute([1, 2, 0])[idx[0]:idx[0] + I.shape[0],
                                   idx[1]:idx[1] + I.shape[1], :]
        PRED += LINEAR
        for c in range(3):
            LINEAR[:, :, :, c] *= medians_whole[c]
            PRED[:, :, :, c] *= medians_whole[c]

        print('\nSharpening / threshdolding / masking')
        # unsharp mask
        blurred = torch.zeros_like(PRED)
        for c in range(PRED.shape[3]):
            blurred[..., c] = gaussian_blur_3d(PRED[..., c], [unsharp_sigma, unsharp_sigma, unsharp_sigma], device)
        PRED += unsharp_amount * (PRED - blurred)
        # mask
        M = (PRED.mean(dim=-1) > threshold_fg)
        for c in range(I.shape[3]):
            PRED[..., c] *= M

    return LINEAR, PRED, affnew

###########

def gaussian_blur_3d(input, stds, device, dtype=torch.float):
    blurred = input[None, None, :, :, :]
    if stds[0] > 0:
        kx = make_gaussian_kernel(stds[0], device=device, dtype=dtype)
        blurred = conv3d(blurred, kx[None, None, :, None, None], stride=1, padding=(len(kx) // 2, 0, 0))
    if stds[1] > 0:
        ky = make_gaussian_kernel(stds[1], device=device, dtype=dtype)
        blurred = conv3d(blurred, ky[None, None, None, :, None], stride=1, padding=(0, len(ky) // 2, 0))
    if stds[2] > 0:
        kz = make_gaussian_kernel(stds[2], device=device, dtype=dtype)
        blurred = conv3d(blurred, kz[None, None, None, None, :], stride=1, padding=(0, 0, len(kz) // 2))
    return torch.squeeze(blurred)

###########

def make_gaussian_kernel(sigma, device, dtype):
    sl = int(np.ceil(3 * sigma))
    ts = torch.linspace(-sl, sl, 2 * sl + 1, dtype=dtype, device=device)
    gauss = torch.exp((-(ts / sigma) ** 2 / 2))
    kernel = gauss / gauss.sum()
    return kernel
