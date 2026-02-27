import os
import numpy as np
from cv2 import imread, resize, distanceTransform
from cv2 import INTER_AREA, DIST_L2
from scipy.io import loadmat
from scipy.ndimage import binary_dilation
import nibabel as nib
from skimage.exposure.exposure import cumulative_distribution
import torch
from scipy.interpolate import RegularGridInterpolator as rgi
import csv


# Read MRI scan
def MRIread(filename, dtype=None, im_only=False):

    assert filename.endswith(
        ('.nii', '.nii.gz', '.mgz')), 'Unknown data file: %s' % filename

    x = nib.load(filename)
    volume = x.get_fdata()
    aff = x.affine

    if dtype is not None:
        volume = volume.astype(dtype=dtype)

    if im_only:
        return volume
    else:
        return volume, aff

# Write MRI scan
def MRIwrite(volume, aff, filename, dtype=None):
    if dtype is not None:
        volume = volume.astype(dtype=dtype)
    if aff is None:
        aff = np.eye(4)
    header = nib.Nifti1Header()
    nifty = nib.Nifti1Image(volume, aff, header)
    nib.save(nifty, filename)

# Crop label volume
def cropLabelVol(V, margin=10, threshold=0):

    # Make sure it's 3D
    margin = np.array(margin)
    if len(margin.shape) < 2:
        margin = [margin, margin, margin]

    if len(V.shape) < 2:
        V = V[..., np.newaxis]
    if len(V.shape) < 3:
        V = V[..., np.newaxis]

    # Now
    idx = np.where(V > threshold)
    i1 = np.max([0, np.min(idx[0]) - margin[0]]).astype('int')
    j1 = np.max([0, np.min(idx[1]) - margin[1]]).astype('int')
    k1 = np.max([0, np.min(idx[2]) - margin[2]]).astype('int')
    i2 = np.min([V.shape[0], np.max(idx[0]) + margin[0] + 1]).astype('int')
    j2 = np.min([V.shape[1], np.max(idx[1]) + margin[1] + 1]).astype('int')
    k2 = np.min([V.shape[2], np.max(idx[2]) + margin[2] + 1]).astype('int')

    cropping = [i1, j1, k1, i2, j2, k2]
    cropped = V[i1:i2, j1:j2, k1:k2]

    return cropped, cropping

# Convert vox to ras coordinates
def vox2ras(vox, vox2ras):
    vox2 = np.concatenate([vox, np.ones(shape=[1, vox.shape[1]])], axis=0)
    ras = np.matmul(vox2ras, vox2)[:-1, :]
    return ras

# Convert ras to vox coordinates
def ras2vox(ras, vox2ras):
    ras2 = np.concatenate([ras, np.ones(shape=[1, ras.shape[1]])], axis=0)
    vox = np.matmul(np.linalg.inv(vox2ras), ras2)[:-1, :]
    return vox

# Apply cropping
def applyCropping(V, cropping):
    i1 = cropping[0]
    j1 = cropping[1]
    k1 = cropping[2]
    i2 = cropping[3]
    j2 = cropping[4]
    k2 = cropping[5]

    if len(V.shape) > 2:
        Vcropped = V[i1:i2, j1:j2, k1:k2, ...]
    else:
        Vcropped = V[i1:i2, j1:j2]

    return Vcropped

# Equalize a grayscale image
def equalize_with_mask(image, mask, nbins=256):
    cdf, bin_centers = cumulative_distribution(image[mask], nbins)
    out = np.interp(image.flat, bin_centers, cdf).reshape(image.shape)
    return (255 * out).astype('uint8')

# Read images and masks from disks, and crop to an equal size so we can easily stack
def read_images_and_masks(d_i, d_s, reverse_ap=False, crop_margin=5, ndilations=5, equalize_images=False):
    print("Extracting slices from photographs")  # try tons of extensions
    Nphotos = len(d_i)
    images = []
    masks = []
    all_croppings = []
    max_x_size = 0
    max_y_size = 0
    total_slice_count = 0
    for n in np.arange(Nphotos):
        X = np.flip(imread(d_i[n]), axis=-1)  # convert to RGB

        if d_s[n][-3:] == "mat":
            Y = loadmat(d_s[n])["LABELS"]
        else:
            Y = np.load(d_s[n])
        print(
            f"Photo {n + 1} has {len(np.unique(Y)) - 1} slices (CCs)"
        )  # Eugenio added -1 to account for zero
        total_slice_count += len(np.unique(Y)) - 1

        for l in 1 + np.arange(np.max(Y)):
            mask, cropping = cropLabelVol(Y == l, crop_margin)
            all_croppings.append(cropping)
            cropping[2] = 0
            cropping[5] = 3
            image = np.squeeze(applyCropping(X, cropping))
            mask = np.squeeze(binary_dilation(mask, iterations=ndilations))
            if equalize_images:
                if len(image.shape)==2:
                    image = equalize_with_mask(image, mask)
                else:
                    for channel in range(3):
                        image[:,:,channel] = equalize_with_mask(image[:,:,channel], mask)
            images.append(image)
            masks.append(mask)
            max_x_size = np.max([max_x_size, mask.shape[0]])
            max_y_size = np.max([max_y_size, mask.shape[1]])

    print(f"Found {total_slice_count} slices in {Nphotos} photos")

    print('Padding the photos to equal size')

    def pad(X, siz):
        if len(X.shape) <= 2:
            X = X[..., None]
        Y = np.zeros([*siz, X.shape[2]])
        idx1 = np.ceil(0.5 * (np.array(siz) - X.shape[:-1])).astype("int")
        idx2 = (idx1 + X.shape[:-1]).astype("int")
        Y[idx1[0]:idx2[0], idx1[1]:idx2[1], :] = X
        return np.squeeze(Y)

    siz = np.round(1.5 * np.array([max_x_size, max_y_size])).astype("int")
    for i in range(total_slice_count):
        images[i] = pad(images[i], siz)
        masks[i] = pad(masks[i], siz)

    if reverse_ap:
        images = images[::-1]
        masks = masks[::-1]

    return images, masks

# Get list of files with photos and masks, in the correct order
def estimate_thicknesses_from_weights(weight_file, areas, av_thickness, max_thickness, reverse_ap=False):
    with open(weight_file, 'r') as file:
        csvreader = csv.reader(file)
        for row in csvreader:
            aux = row.copy()
        weights = np.array(aux, dtype=np.float32)
    if len(weights) != (len(areas) + 1):
        raise Exception('Number of weights is not equal to the number of slices plus one')
    if reverse_ap:
        weights = weights[::-1]
    # thicknesses: trapezoid model for all slabs, except for first and last (ellipsoids)
    thicknesses = np.zeros_like(weights)
    thicknesses[1:-1] = 2.0 * weights[1:-1] / (areas[0:-1] +  areas[1:])
    thicknesses[0] = weights[0] * 3.0 / (2  * areas[0])
    thicknesses[-1] = weights[-1] * 3.0 / (2 * areas[-1])
    thicknesses *= (av_thickness / thicknesses.mean())

    if max_thickness<=0: # not provided
        if np.any(thicknesses > (2.0 * av_thickness)):
            print('*** ERROR ***')
            print('After thickness estimation, one slice is more than twice as thick as the (provided) average')
            print('The (provided) weights and estimated thicknesses are (asterisks mark offenders): ')
            for i in range(len(weights)):
                tline = 'Slab ' + str(i+1) + ', weight = ' + str(weights[i]) + ', thickness = ' + str(thicknesses[i])
                if thicknesses[i] > (2.0 * av_thickness):
                    tline = ' * ' + tline
                print(tline)
            print('This could be an error in the weights (please double check),')
            print('or an overestimation by the geometric model.')
            print('In the latter case, please rerun with the option --thickness_cap set to')
            print('what you believe the maximum thickness should be')
            raise Exception('Exiting...')
    else: # thickness cap provided
        print('Applying thickness cap: ' + str(max_thickness) + ' mm')
        thicknesses[thicknesses > max_thickness] = max_thickness
        thicknesses *= (av_thickness / thicknesses.mean())

    return thicknesses


def make_distance_transforms(mask_list, pixel_size):
    distances = []
    for mask in mask_list:
        D = distanceTransform(np.uint8(mask > 0), DIST_L2, 5).astype(np.float32) - 0.5
        Din = 0.5 - distanceTransform(np.uint8(mask == 0), DIST_L2, 5).astype(np.float32)
        D[mask==0] = Din[mask==0]
        distances.append(pixel_size * D)
    return distances

# Resample images and masks to a specified resolution
def resample_images_and_masks(images, masks, photo_res, output_res, stretch_lr=1.0):
    print("Resampling to the output resolution: " + str(output_res) + " mm")
    Nslices = len(images)
    images_resized = []
    masks_resized = []
    for n in np.arange(Nslices):
        Isl = resize(
            images[n],
            None,
            fx=photo_res * stretch_lr / output_res,
            fy=photo_res / output_res,
            interpolation=INTER_AREA,
        )
        Msl = resize(
            masks[n].astype("float"),
            None,
            fx=photo_res * stretch_lr / output_res,
            fy=photo_res / output_res,
            interpolation=INTER_AREA,
        )

        images_resized.append(Isl)
        masks_resized.append(Msl)

    return images_resized, masks_resized

# Make resolution pyramid for images and masks
def make_pyramid(images, masks, resolutions, original_res, slice_thickness, reverse_lr):

    print("Making resolution pyramid")
    Is = []  # images
    Ms = []  # masks
    Affs = []  # affine matrices (a.k.a. "voxel-to-RAS" or "vox2ras")
    Nslices = images.shape[2]
    Nscales = len(resolutions)

    # Make affine matrix of the inputStart with output resolution: affine matrix and padding
    aff_orig = np.array(
        [
            [0, -original_res, 0, 0],
            [0, 0, -slice_thickness, 0],
            [-original_res, 0, 0, 0],
            [0, 0, 0, 1],
        ]
    )
    if reverse_lr:
        aff_orig[0, 1] = -aff_orig[0, 1]

    # go over scales
    for s in np.arange(Nscales - 1, -1, -1):

        # let's take care of the affine matrix first
        aff = np.zeros([4, 4])
        aff[0, 1] = -resolutions[s]
        aff[1, 2] = -slice_thickness
        aff[2, 0] = -resolutions[s]
        aff[3, 3] = 1
        if reverse_lr:
            aff[0, 1] = -aff[0, 1]
        aux = np.array([[resolutions[s] / original_res], [resolutions[s] / original_res], [1],])
        aff[0:3, 3] = np.matmul(aff_orig[0:3, 0:3], (0.5 * aux - 0.5))[:, 0]

        # go over slices
        for n in range(Nslices):
            f = (original_res / resolutions[s])
            Isl = resize(np.mean(images[:, :, n, :], axis=-1), None, fx=f, fy=f, interpolation=INTER_AREA)
            Msl = resize(masks[:, :, n], None, fx=f, fy=f, interpolation=INTER_AREA)
            if n == 0:
                im = np.zeros([*Msl.shape, Nslices])
                mask = np.zeros([*Msl.shape, Nslices])
            im[:, :, n] = Isl
            mask[:, :, n] = Msl

        Is.insert(0, im)
        Ms.insert(0, mask)
        Affs.insert(0, aff)

    return Is, Ms, Affs, aff_orig

# Compute hemisphere mask from SynthSeg
def compute_mask_from_synthseg(SEG, hemi):
    lut=np.ones(3000, dtype=np.int32)
    # Kill cerebellu, brainstem, ventricles, CSF
    lut[0] = lut[7:9] = lut[46:48] = lut[4:6] = lut[43:45] = lut[24] = lut[14:17] = 0
    if hemi=='left':
        lut[35:70] = lut[2000:] = 0
    if hemi=='right':
        lut[:35] = lut[1000:2000] = 0
    MASK = lut[SEG.round().astype(np.int32)]
    return MASK

# Center affine matrix in origin of RAS coordinate system using centroid of data matrix
def center_affine_matrix_in_origin(X, aff):
    idx = np.where(X > 0)
    cog_vox = np.array(
        [[np.mean(idx[0])], [np.mean(idx[1])], [np.mean(idx[2])]]
    )
    cog_ras = vox2ras(cog_vox, aff)
    aff_centered = np.copy(aff)
    aff_centered[:-1, -1] = aff_centered[:-1, -1] - np.squeeze(cog_ras)
    return aff_centered, cog_ras

# Nearst negithbor or trilinear 3D interpolation with pytorch
def fast_3D_interp_torch(X, II, JJ, KK, mode):
    device = X.device
    if mode=='nearest':
        IIr = torch.round(II).long()
        JJr = torch.round(JJ).long()
        KKr = torch.round(KK).long()
        IIr[IIr < 0] = 0
        JJr[JJr < 0] = 0
        KKr[KKr < 0] = 0
        IIr[IIr > (X.shape[0] - 1)] = (X.shape[0] - 1)
        JJr[JJr > (X.shape[1] - 1)] = (X.shape[1] - 1)
        KKr[KKr > (X.shape[2] - 1)] = (X.shape[2] - 1)
        if len(X.shape)==3:
            X = X[..., None]
        Y = torch.zeros([*II.shape, X.shape[3]], device=device)
        for channel in range(X.shape[3]):
            aux = X[:, :, :, channel]
            Y[:,:,:,channel] = aux[IIr, JJr, KKr]
        if Y.shape[3] == 1:
            Y = Y[:, :, :, 0]

    elif mode=='linear':
        ok = (II>0) & (JJ>0) & (KK>0) & (II<=X.shape[0]-1) & (JJ<=X.shape[1]-1) & (KK<=X.shape[2]-1)
        IIv = II[ok]
        JJv = JJ[ok]
        KKv = KK[ok]

        fx = torch.floor(IIv).long()
        cx = fx + 1
        cx[cx > (X.shape[0] - 1)] = (X.shape[0] - 1)
        wcx = IIv - fx
        wfx = 1 - wcx

        fy = torch.floor(JJv).long()
        cy = fy + 1
        cy[cy > (X.shape[1] - 1)] = (X.shape[1] - 1)
        wcy = JJv - fy
        wfy = 1 - wcy

        fz = torch.floor(KKv).long()
        cz = fz + 1
        cz[cz > (X.shape[2] - 1)] = (X.shape[2] - 1)
        wcz = KKv - fz
        wfz = 1 - wcz

        if len(X.shape)==3:
            X = X[..., None]

        Y = torch.zeros([*II.shape, X.shape[3]], device=device)
        for channel in range(X.shape[3]):
            Xc = X[:, :, :, channel]

            c000 = Xc[fx, fy, fz]
            c100 = Xc[cx, fy, fz]
            c010 = Xc[fx, cy, fz]
            c110 = Xc[cx, cy, fz]
            c001 = Xc[fx, fy, cz]
            c101 = Xc[cx, fy, cz]
            c011 = Xc[fx, cy, cz]
            c111 = Xc[cx, cy, cz]

            c00 = c000 * wfx + c100 * wcx
            c01 = c001 * wfx + c101 * wcx
            c10 = c010 * wfx + c110 * wcx
            c11 = c011 * wfx + c111 * wcx

            c0 = c00 * wfy + c10 * wcy
            c1 = c01 * wfy + c11 * wcy

            c = c0 * wfz + c1 * wcz

            Yc = torch.zeros(II.shape, device=device)
            Yc[ok] = c.float()
            Y[:,:,:,channel] = Yc

        if Y.shape[3]==1:
            Y = Y[:,:,:,0]

    else:
        raise Exception('mode must be linear or nearest')

    return Y

# Computes deformation as a field of RAS coordinates, in photo space and at 1mm isotropic
def computeRAS(grids_new_mri_nonlin, REFshape, REFaff, photo_aff, iso_size, iso_aff, device):
    # First in native space of the photo recon
    IJKmri = []
    for c in range(3):
        IJKmri.append(0.5 * ((grids_new_mri_nonlin[:, :, :, 2 - c] + 1) * (REFshape[c] - 1)))
    IJKmri = np.stack(IJKmri, axis=-1)
    RAS = np.zeros_like(IJKmri)
    for c in range(3):
        RAS[..., c] = REFaff[c, 0] * IJKmri[..., 0] + REFaff[c, 1] * IJKmri[..., 1] + REFaff[c, 2] * IJKmri[
            ..., 2] + REFaff[c, 3]
    # now in the 1x1x1 space of the linear / ML interpolations
    v2v = torch.tensor(np.linalg.inv(photo_aff) @ iso_aff, device=device, dtype=torch.float32)
    coords = torch.stack(torch.meshgrid(
        torch.arange(iso_size[0], device=device),
        torch.arange(iso_size[1], device=device),
        torch.arange(iso_size[2], device=device),
        indexing='ij'), dim=-1).reshape(-1, 3)  # shape (N, 3)
    coords_trans = (v2v @ torch.cat([coords, torch.ones(coords.shape[0], 1, device=device)], dim=-1).T).T[:,
                   :3]  # (N,3)
    grid = 2.0 * coords_trans / (torch.tensor(RAS.shape[:3], device=device) - 1) - 1.0  # scale to [-1,1]
    grid = grid.view(1, iso_size[0], iso_size[1], iso_size[2], 3)  # shape (1, X, Y, Z, 3)
    RASres = F.grid_sample(torch.tensor(RAS, device=device, dtype=torch.float32).permute([3, 0, 1, 2]).unsqueeze(0),
                           grid.flip([-1]), mode='bilinear', padding_mode='border', align_corners=True)
    RASres = RASres[0].permute(1, 2, 3, 0).detach().cpu().numpy()
    return RAS, RASres

# deform an image with a field of RAS coordinates
def deform(I, aff, RAS, device, mode='linear'):
    aff_inv = np.linalg.inv(aff)
    vox = np.zeros_like(RAS)
    for c in range(3):
        vox[..., c] = aff_inv[c, 3]
        for cc in range(3):
            vox[..., c] += aff_inv[c, cc] * RAS[..., cc]
    vt = torch.tensor(vox).to(device)
    D = fast_3D_interp_torch(torch.tensor(I).to(device), vt[..., 0], vt[..., 1], vt[..., 2], mode).cpu().detach().numpy()
    return D, vox

# Deforms volumes and surfaces in FS directory
def deform_FS_derivatives(deform_recon_dir, hemisphere, RAS, RASaff, output_directory, device):

    print('  Deforming  aparc+aseg')
    aaseg, aaseg_aff = MRIread(deform_recon_dir + '/mri/aparc+aseg.mgz')
    mask = compute_mask_from_synthseg(aaseg, hemisphere)
    aaseg *= mask
    AASEGresampled, vox = deform(aaseg, aaseg_aff, RAS, device, mode='nearest')
    MRIwrite(AASEGresampled, RASaff, output_directory + '/aparc+aseg.deformed.mgz')

    print('  Inverting field (approximate) to deform surfaces')
    # We use kernel regression with a tiny kernel (one voxel, weights of trilinear interpolation)
    X = np.zeros([*aaseg.shape, 3])
    XdenN = np.zeros(aaseg.shape)
    II = vox[...,0]; JJ = vox[..., 1]; KK = vox[..., 2];
    ok = (II > 0) & (JJ > 0) & (KK > 0) & (II <= X.shape[0] - 1) & (JJ <= X.shape[1] - 1) & (KK <= X.shape[2] - 1)
    IIv = II[ok]; JJv = JJ[ok]; KKv = KK[ok]
    grid = torch.stack(torch.meshgrid([torch.arange(RAS.shape[s]) for s in range(3)]), dim=-1).cpu().detach().numpy()
    rs = []
    for c in range(3):
        rs.append((RASaff[c, 0] * grid[..., 0] + RASaff[c, 1] * grid[..., 1] + RASaff[c, 2] * grid[..., 2] + RASaff[c, 3])[ok])
    rs = np.stack(rs, axis=1)

    fx = np.floor(IIv).astype(np.int32); cx = fx + 1; cx[cx > (X.shape[0] - 1)] = (X.shape[0] - 1)
    fy = np.floor(JJv).astype(np.int32); cy = fy + 1; cy[cy > (X.shape[1] - 1)] = (X.shape[1] - 1)
    fz = np.floor(KKv).astype(np.int32); cz = fz + 1; cz[cz > (X.shape[2] - 1)] = (X.shape[2] - 1)
    wcx = IIv - fx; wfx = 1 - wcx;
    wcy = JJv - fy; wfy = 1 - wcy;
    wcz = KKv - fz; wfz = 1 - wcz;
    fff = (wfx * wfy * wfz)
    ffc = (wfx * wfy * wcz)
    fcf = (wfx * wcy * wfz)
    fcc = (wfx * wcy * wcz)
    cff = (wcx * wfy * wfz)
    cfc = (wcx * wfy * wcz)
    ccf = (wcx * wcy * wfz)
    ccc = (wcx * wcy * wcz)
    for idx, (fx_,fy_,fz_,cx_,cy_,cz_,fff_,ffc_,fcf_,fcc_,cff_,cfc_,ccf_,ccc_,rs_) in enumerate(zip(fx,fy,fz,cx,cy,cz,fff,ffc,fcf,fcc,cff,cfc,ccf,ccc,rs)):
        if (idx%1000000)==0:
            print('Voxel ' + str(idx) + ' of ' + str(len(fx)))
        XdenN[fx_, fy_, fz_] += fff_
        XdenN[fx_, fy_, cz_] += ffc_
        XdenN[fx_, cy_, fz_] += fcf_
        XdenN[fx_, cy_, cz_] += fcc_
        XdenN[cx_, fy_, fz_] += cff_
        XdenN[cx_, fy_, cz_] += cfc_
        XdenN[cx_, cy_, fz_] += ccf_
        XdenN[cx_, cy_, cz_] += ccc_
        X[fx_, fy_, fz_, :] += (fff_ * rs_)
        X[fx_, fy_, cz_, :] += (ffc_ * rs_)
        X[fx_, cy_, fz_, :] += (fcf_ * rs_)
        X[fx_, cy_, cz_, :] += (fcc_ * rs_)
        X[cx_, fy_, fz_, :] += (cff_ * rs_)
        X[cx_, fy_, cz_, :] += (cfc_ * rs_)
        X[cx_, cy_, fz_, :] += (ccf_ * rs_)
        X[cx_, cy_, cz_, :] += (ccc_ * rs_)

    filter = torch.zeros([3,3,3]).to(device)
    filter[1,1,:]=0.05; filter[:,1,1]=0.05; filter[1,:,1]=0.05; filter[1,1,1]=0.7
    X = torch.tensor(X,dtype=filter.dtype).to(device)
    XdenN = torch.tensor(XdenN,dtype=filter.dtype).to(device)
    stride = (1, 1, 1); padding = (1,1,1)
    XdenN = torch.nn.functional.conv3d(XdenN[None, None, ...], filter[None, None], stride=stride, padding=padding)[0, 0, ...]
    X = torch.nn.functional.conv3d(X[None,...].permute([4,0,1,2,3]), filter[None, None], stride=stride, padding=padding).squeeze().permute([1,2,3,0])
    X /= (1e-6 + XdenN[...,None])
    MRIwrite(X.cpu().detach().numpy(), aaseg_aff, output_directory + '/ras.inv.nii.gz')

    print('  Deforming surfaces')
    interp = rgi((np.arange(X.shape[0]), np.arange(X.shape[1]), np.arange(X.shape[2])), X.cpu().detach().numpy())
    file_list = ['lh.white', 'lh.pial', 'rh.white', 'rh.pial']
    if hemisphere == 'right':
        file_list = ['rh.white', 'rh.pial']
    if hemisphere == 'left':
        file_list = ['lh.white', 'lh.pial']
    for surf_file in file_list:
        coords, faces, meta = nib.freesurfer.io.read_geometry(deform_recon_dir + '/surf/' + surf_file, read_metadata=True)
        coords += meta['cras'][None, :]
        meta['cras'][:] = 0
        vox_norm = np.linalg.inv(aaseg_aff) @ np.concatenate([coords, np.ones([coords.shape[0], 1])],  axis=1).transpose()
        coords_deformed = interp(vox_norm[:-1, :].T)
        nib.freesurfer.io.write_geometry(output_directory + '/' + surf_file + '.deformed', coords_deformed, faces, volume_info=meta)


# Integrate a bunch of 2D SVFs into deformation fields
def svf_integration(svf, n_steps=6):
    # Rearrange to (N, 2, H, W) for grid_sample
    v = svf.permute(3, 0, 1, 2).contiguous()
    N, C, H, W = v.shape
    v = v / (2 ** n_steps)
    r = torch.linspace(-1, 1, H, device=v.device)
    c = torch.linspace(-1, 1, W, device=v.device)
    grid_r, grid_c = torch.meshgrid(r, c, indexing='ij')
    base_grid = torch.stack((grid_r, grid_c), dim=-1)  # (H, W, 2)
    base_grid = base_grid.unsqueeze(0).repeat(N, 1, 1, 1)  # (N, H, W, 2)

    # Convert displacement from pixel space to normalized grid space
    v_norm = torch.zeros_like(v)
    v_norm[:, 0] = v[:, 0] * (2.0 / (H - 1))
    v_norm[:, 1] = v[:, 1] * (2.0 / (W - 1))

    # Initialize deformation field φ = Id + v
    phi = base_grid + v_norm.permute(0, 2, 3, 1)

    # Scaling and squaring
    for _ in range(n_steps):
        # Sample phi at phi locations (compose deformation)
        phi = F.grid_sample(
            phi.permute(0, 3, 1, 2),
            phi.flip([3]),
            mode='bilinear',
            padding_mode='border',
            align_corners=True
        ).permute(0, 2, 3, 1)

    # Convert back to displacement field
    disp = phi - base_grid
    disp_pix = torch.zeros_like(disp)
    disp_pix[..., 0] = disp[..., 0] * (H - 1) / 2.0
    disp_pix[..., 1] = disp[..., 1] * (W - 1) / 2.0

    # Rearrange back to (2, H, W, N)
    disp_pix = disp_pix.permute(3, 1, 2, 0).contiguous()

    return disp_pix
