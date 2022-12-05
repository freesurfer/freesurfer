# Imports
import argparse
import glob
import os
import sys

import cv2
import nibabel as nib
import numpy as np
import scipy.io
import scipy.ndimage
import torch
import torch.nn as nn
import torch.nn.functional as nnf
from torch.optim import Optimizer

from copy import deepcopy
from functools import reduce

# ================================================================================================
#                                         Main Entrypoint
# ================================================================================================


def main():

    torch.set_default_dtype(torch.float64)
    seed_all(0)

    ########################################################
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Code for 3D photo reconstruction (Tregidgo, et al., MICCAI 2020)"
    )

    parser.add_argument(
        "--input_photo_dir",
        type=str,
        nargs="*",
        help="Directory with input photos (required)",
        required=True,
    )
    parser.add_argument(
        "--input_segmentation_dir",
        type=str,
        nargs="*",
        help="Directory with input slab masks / segmentations (required)",
        required=True,
    )

    parser.add_argument(
        "--ref_mask", type=str, help="Reference binary mask", default=None
    )
    parser.add_argument(
        "--ref_surface", type=str, help="Reference surface file", default=None
    )
    parser.add_argument(
        "--ref_soft_mask", type=str, help="Reference soft mask", default=None
    )

    parser.add_argument(
        "--mesh_reorient_with_indices",
        type=str,
        help="Vertex indices of frontal pole, occipital pole, and top of central sulcus, separated with commas, for mesh alignment",
        default=None,
    )

    parser.add_argument(
        "--photos_of_posterior_side",
        dest="posterior_side",
        action="store_true",
        help="Use when photos are taken of posterior side of slabs (default is anterior side)",
    )
    parser.set_defaults(posterior_side=False)

    parser.add_argument(
        "--order_posterior_to_anterior",
        dest="posterior_to_anterior",
        action="store_true",
        help="Use when photos are ordered from posterior to anterior (default is anterior to posterior)",
    )
    parser.set_defaults(posterior_to_anterior=False)

    parser.add_argument(
        "--allow_z_stretch",
        dest="allow_z_stretch",
        action="store_true",
        help="Use to adjust the slice thickness to best match the reference."
        + " You should probably *never* use this with soft references (ref_soft_mask)",
    )
    parser.set_defaults(allow_z_stretch=False)

    parser.add_argument(
        "--rigid_only_for_photos",
        dest="rigid_only_for_photos",
        action="store_true",
        help="Switch on if you want photos to deform only rigidly (not affine)",
    )
    parser.set_defaults(rigid_only_for_photos=False)

    parser.add_argument(
        "--slice_thickness", type=float, help="Slice thickness in mm", required=True
    )
    parser.add_argument(
        "--photo_resolution",
        type=float,
        help="Resolution of the photos in mm",
        required=True,
    )

    parser.add_argument(
        "--output_directory",
        type=str,
        help="Output directory with reconstructed photo volume and reference",
        required=True,
    )

    parser.add_argument("--gpu", type=int, help="Index of GPU to use", default=None)

    options = parser.parse_args()

    ########################################################
    # Set the GPU if needed
    if options.gpu is None:
        print("Using the CPU")
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
        device = torch.device("cpu")
    else:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(options.gpu)
        if torch.cuda.is_available():
            print("Using GPU device " + str(options.gpu))
            device = torch.device("cuda:0")
        else:
            print(
                "Tried to use GPU device "
                + str(options.gpu)
                + " but failed; using CPU instead"
            )
            device = torch.device("cpu")

    ########################################################
    # Input data

    # First, make sure that you specify one and only one reference!
    n_refs = (
        (options.ref_mask is not None)
        + (options.ref_soft_mask is not None)
        + (options.ref_surface is not None)
    )
    if n_refs != 1:
        raise Exception(
            "You should provide 1 and only 1 reference: binary mask, soft mask, or surface"
        )

    if options.ref_mask is not None:
        input_reference = options.ref_mask
        ref_type = "mask"
    elif options.ref_soft_mask is not None:
        input_reference = options.ref_soft_mask
        ref_type = "soft_mask"
    else:
        input_reference = options.ref_surface
        ref_type = "surface"

    if not os.path.exists(input_reference):
        sys.exit(f"DNE: {input_reference}")

    # Some directory names have spaces... oh well
    input_photo_dir = options.input_photo_dir[0]
    for i in range(len(options.input_photo_dir) - 1):
        input_photo_dir = input_photo_dir + " " + options.input_photo_dir[i + 1]
    input_photo_dir = input_photo_dir + "/"

    input_segmentation_dir = options.input_segmentation_dir[0]
    for i in range(len(options.input_segmentation_dir) - 1):
        input_segmentation_dir = (
            input_segmentation_dir + " " + options.input_segmentation_dir[i + 1]
        )
    input_segmentation_dir = input_segmentation_dir + "/"

    reverse_lr = options.posterior_side
    reverse_ap = options.posterior_to_anterior
    slice_thickness = options.slice_thickness
    photo_res = options.photo_resolution
    allow_z_stretch = options.allow_z_stretch

    # Outputs
    output_directory = options.output_directory + "/"
    output_photo_recon = output_directory + "photo_recon.mgz"
    if ref_type == "surface":
        output_registered_reference = output_directory + "registered_reference.surf"
    else:
        output_registered_reference = None

    if os.path.isdir(output_directory) is False:
        os.makedirs(output_directory, exist_ok=True)

    # Training
    alternate_optimization = True

    ########################################################

    # Constants
    RESOLUTIONS = [4, 2, 1, 0.5]
    STEPS = [300, 300, 300, 300]
    if RESOLUTIONS[-1] < photo_res:
        RESOLUTIONS[-1] = photo_res

    LR = 1.0
    TOL = 1e-6

    if ref_type == "mask":
        K_DICE_MRI = 0.95
        K_DICE_SLICES = 0.025
        K_NCC_SLICES = 0.025
        K_REGULARIZER = 0.00025
        K_NCC_INTERMODALITY = None
        K_SURFACE_TERM = None
    elif ref_type == "soft_mask":
        K_DICE_MRI = 0.8
        K_DICE_SLICES = 0.1
        K_NCC_SLICES = 0.1
        K_REGULARIZER = 0.1
        K_NCC_INTERMODALITY = None
        K_SURFACE_TERM = None
    else:  # surface
        K_SURFACE_TERM = 1.0
        K_DICE_SLICES = 0.025
        K_NCC_SLICES = 0.025
        K_REGULARIZER = 0.0025
        K_DICE_MRI = 0.95
        K_NCC_INTERMODALITY = None


    ########################################################

    print("Extracting slices from photographs")
    d_i = glob.glob(input_photo_dir + "/*.jpg")
    if len(d_i) == 0:
        d_i = glob.glob(input_photo_dir + "/*.tif")
    if len(d_i) == 0:
        d_i = glob.glob(input_photo_dir + "/*.tiff")
    if len(d_i) == 0:
        d_i = glob.glob(input_photo_dir + "/*.JPG")
    if len(d_i) == 0:
        d_i = glob.glob(input_photo_dir + "/*.png")
    d_i = sorted(d_i)

    d_s = glob.glob(input_segmentation_dir + "/*.mat")
    if len(d_s) == 0:
        d_s = glob.glob(input_segmentation_dir + "/*.npy")
    d_s = sorted(d_s)

    if reverse_ap:
        d_s = d_s[::-1]
        d_i = d_i[::-1]

    Nphotos = len(d_i)

    Iorig = []
    Morig = []

    all_croppings = []
    total_slice_count = 0
    for n in np.arange(Nphotos):
        X = np.flip(cv2.imread(d_i[n]), axis=-1)  # convert to RGB

        if d_s[n][-3:] == "mat":
            Y = scipy.io.loadmat(d_s[n])["LABELS"]
        else:
            Y = np.load(d_s[n])
        print(
            f"Photo {n + 1} has {len(np.unique(Y))-1} slices (CCs)"
        )  # Eugenio added -1 to account for zero
        total_slice_count += len(np.unique(Y)) - 1

        for l in 1 + np.arange(np.max(Y)):
            mask, cropping = cropLabelVol(Y == l, np.round(5 / photo_res))
            all_croppings.append(cropping)
            cropping[2] = 0
            cropping[5] = 3
            image = applyCropping(X, cropping)
            image = image * mask
            Iorig.append(image)
            Morig.append(mask)

    print(f"Found {total_slice_count} slices in {Nphotos} photos")

    ########################################################

    print(
        "Resampling to highest target resolution: " + str(RESOLUTIONS[-1]) + " mm"
    )

    Nslices = len(Iorig)
    Nscales = len(RESOLUTIONS)
    I = []
    M = []
    for n in np.arange(Nslices):
        Isl = cv2.resize(
            Iorig[n],
            None,
            fx=photo_res / RESOLUTIONS[-1],
            fy=photo_res / RESOLUTIONS[-1],
            interpolation=cv2.INTER_AREA,
        )
        Msl = cv2.resize(
            Morig[n].astype("float"),
            None,
            fx=photo_res / RESOLUTIONS[-1],
            fy=photo_res / RESOLUTIONS[-1],
            interpolation=cv2.INTER_AREA,
        )
        Isl[Msl == 0] = 0
        I.append(Isl)
        M.append(Msl)

    ########################################################

    print("Coarse alignment and padding")

    PAD_AP = 3  # Padding in A-P axis ensures that mask remains 100% in FOV
    sizes = np.zeros([Nslices, 2])
    for n in np.arange(Nslices):
        sizes[n, :] = M[n].shape
    siz = np.round(1.5 * np.max(sizes, axis=0)).astype("int")
    Is = []
    Ms = []
    Affs = []

    aff = np.array(
        [
            [0, -RESOLUTIONS[-1], 0, 0],
            [0, 0, -slice_thickness, 0],
            [-RESOLUTIONS[-1], 0, 0, 0],
            [0, 0, 0, 1],
        ]
    )
    if reverse_lr:
        aff[0, 1] = -aff[0, 1]
    im = np.zeros([*siz, Nslices + 2 * PAD_AP, 3])
    mask = np.zeros([*siz, Nslices + 2 * PAD_AP])
    all_paddings = []
    for n in np.arange(Nslices):
        idx1 = np.ceil(0.5 * (np.array(siz) - sizes[n, :])).astype("int")
        idx2 = (idx1 + sizes[n, :]).astype("int")
        im[idx1[0] : idx2[0], idx1[1] : idx2[1], n + PAD_AP, :] = I[n]
        mask[idx1[0] : idx2[0], idx1[1] : idx2[1], n + PAD_AP] = M[n]
        all_paddings.append(idx1[0:2] - all_croppings[n][0:2])

    Is.append(im)
    Ms.append(mask)
    Affs.append(aff)
    Nslices = Nslices + 2 * PAD_AP

    ########################################################

    print("Building resolution pyramid")

    for s in np.arange(Nscales - 2, -1, -1):
        for n in np.arange(Nslices):
            Isl = cv2.resize(
                Is[-1][:, :, n, :],
                None,
                fx=RESOLUTIONS[-1] / RESOLUTIONS[s],
                fy=RESOLUTIONS[-1] / RESOLUTIONS[s],
                interpolation=cv2.INTER_AREA,
            )
            Msl = cv2.resize(
                Ms[-1][:, :, n],
                None,
                fx=RESOLUTIONS[-1] / RESOLUTIONS[s],
                fy=RESOLUTIONS[-1] / RESOLUTIONS[s],
                interpolation=cv2.INTER_AREA,
            )
            if n == 0:
                im = np.zeros([*Msl.shape, Nslices, 3])
                mask = np.zeros([*Msl.shape, Nslices])

            im[:, :, n, :] = Isl
            mask[:, :, n] = Msl

        aff = np.zeros([4, 4])
        aff[0, 1] = -RESOLUTIONS[s]
        aff[1, 2] = -slice_thickness
        aff[2, 0] = -RESOLUTIONS[s]
        aff[3, 3] = 1
        if reverse_lr:
            aff[0, 1] = -aff[0, 1]
        aux = np.array(
            [
                [RESOLUTIONS[s] / RESOLUTIONS[-1]],
                [RESOLUTIONS[s] / RESOLUTIONS[-1]],
                [1],
            ]
        )
        aff[0:3, 3] = np.matmul(Affs[-1][0:3, 0:3], (0.5 * aux - 0.5))[:, 0]

        Is.insert(0, im)
        Ms.insert(0, mask)
        Affs.insert(0, aff)

    # Switch from intensities to gradients, if we are using surfaces
    Is_copy = np.copy(Is[-1])
    if ref_type == "surface":
        for s in range(Nscales):
            erode_its = np.ceil(1.0 / RESOLUTIONS[s]).astype("int")
            for z in range(Nslices):
                M_ERODED = scipy.ndimage.binary_erosion(
                    Ms[s][:, :, z] > 0.5, iterations=erode_its
                )
                for c in range(3):
                    Is[s][:, :, z, c] = grad2d(Is[s][:, :, z, c]) * M_ERODED

    ########################################################

    print("Reading and preprocessing reference")

    if ref_type != "surface":
        REF, REFaff = MRIread(input_reference)
        if np.isnan(REF).any():
            print("There are NaNs here")
            REF[np.isnan(REF)] = 0
        REF = np.squeeze(REF)

        # Eugenio: added padding here; usefull when using hard masks from mris_fill
        pad = 10
        REF_padded = np.zeros(np.array(REF.shape) + 2 * pad)
        REF_padded[pad:-pad, pad:-pad, pad:-pad] = REF
        REFaff_padded = np.copy(REFaff)
        REFaff_padded[:-1, -1] = REFaff_padded[:-1, -1] - np.matmul(
            REFaff_padded[:-1, :-1], np.array([pad, pad, pad])
        )
        REF = REF_padded
        REFaff = REFaff_padded

        REF_orig = np.copy(REF)
        REFaff_orig = np.copy(REFaff)

        REF = np.squeeze(REF) / np.max(REF)
        if ref_type == "mask" or ref_type == "image":
            REF = (REF > 0).astype("float")

    # Surfaces require quite a bit of extra work
    else:

        input_mesh_converted = output_directory + "/input_mesh.surf"
        input_mesh_reoriented = output_directory + "/input_mesh_reoriented.surf"
        reoriented_mask_vol = output_directory + "/input_mesh_reoriented.filled.mgz"

        fs_home = os.getenv("FREESURFER_HOME")
        if fs_home is None:
            raise Exception(
                "FREESURFER_HOME variable not found; is FreeSurfer sourced?"
            )

        print("Converting reference mesh to FreeSurfer format")
        a = os.system(
            "mris_convert "
            + input_reference
            + " "
            + input_mesh_converted
            + " >/dev/null"
        )
        if a > 0:
            raise Exception("error in mris_convert... is FreeSurfer sourced?")

        print()
        # Read in and fill in missing metadata if needed (eg if STL file)
        P, T, meta = nib.freesurfer.read_geometry(
            input_mesh_converted, read_metadata=True
        )
        if meta["valid"][0] == "0":
            meta["valid"] = "1  # volume info valid"
            meta["filename"] = ""
            meta["volume"] = np.array([256, 256, 256]).astype(int)
            meta["voxelsize"] = np.array([1.0, 1.0, 1.0])
            meta["xras"] = np.array([-1.0, 0.0, 0.0])
            meta["yras"] = np.array([0.0, 0.0, -1.0])
            meta["zras"] = np.array([0.0, 1.0, 0.0])
            meta["cras"] = np.array([0.0, 0.0, 0.0])

        # Apply rotation using provided key vertices, if provided
        # https://towardsdatascience.com/the-definitive-procedure-for-aligning-two-sets-of-3d-points-with-the-kabsch-algorithm-a7ec2126c87e
        if options.mesh_reorient_with_indices is None:
            print("No indices were provided to reorient mesh; just copying over...")
            a = os.system(
                "cp "
                + input_mesh_converted
                + " "
                + input_mesh_reoriented
                + " >/dev/null "
            )
            if a > 0:
                raise Exception("error copying mesh")
        else:
            print("Reorienting mesh with provided vertices")
            idx = np.zeros(3).astype(int)
            aux = options.mesh_reorient_with_indices.split(",")
            for i in range(len(idx)):
                idx[i] = int(aux[i])
            K = P[idx, :]
            K = K - np.mean(K, axis=0)

            if True:  # rough RAS aligment, already demeaned!
                Kref = np.array([[0, 85, -20], [0, -80, -25], [0, -5, 45]]).astype(
                    float
                )
            else:  # precomputed from rh.white
                Kref = np.array(
                    [
                        [5.64194918, 77.57227325, 10.32956219],
                        [1.60726917, -90.65991211, -0.76444769],
                        [3.86025476, -13.81834793, 69.90812683],
                    ]
                )
                Kref = Kref - np.mean(Kref, axis=0)

            H = np.transpose(Kref) @ K
            U, S, Vt = np.linalg.svd(H)
            if np.linalg.det(np.transpose(Vt) @ U) > 0:
                R = np.transpose(Vt) @ np.transpose(U)
            else:
                E = np.eye(3)
                E[2, 2] = -1
                R = np.transpose(Vt) @ (E @ np.transpose(U))

            P = P - np.mean(P, axis=0)
            P = P @ R
            meta["cras"][:] = 0
            nib.freesurfer.write_geometry(
                input_mesh_reoriented, P, T, volume_info=meta
            )

        # Fill in the mesh
        print("Filling in mesh to obtain binary volume")
        a = os.system(
            "mris_fill -r 1 "
            + input_mesh_reoriented
            + " "
            + output_directory
            + "/temp.mgz >/dev/null"
        )
        if a > 0:
            raise Exception("error in mris_fill... is FreeSurfer sourced?")
        # We pad a bit
        [img, aff] = MRIread(output_directory + "/temp.mgz")
        pad = 8
        img2 = np.zeros(np.array(img.shape) + 2 * pad)
        img2[pad:-pad, pad:-pad, pad:-pad] = img
        aff2 = aff
        aff2[:-1, -1] = aff2[:-1, -1] - np.squeeze(
            np.matmul(aff[:-1, :-1], pad * np.ones([3, 1]))
        )
        MRIwrite(img2, aff2, reoriented_mask_vol)
        os.system("rm -rf " + output_directory + "/temp.mgz >/dev/null")

        # Read deformed surface and corresponding reference volume
        REF, REFaff = MRIread(reoriented_mask_vol)
        REF_orig = np.copy(REF)
        REFaff_orig = np.copy(REFaff)
        REF = np.squeeze(REF) / np.max(REF)
        REF = (REF > 0.5).astype("float")

        Pmesh, Tmesh, meta_mesh = nib.freesurfer.read_geometry(
            input_mesh_reoriented, read_metadata=True
        )
        Pmesh += meta_mesh["cras"]  # ** CRUCIAL **
        meta_mesh["cras"][:] = 0

        # And finally, take the gradient of the photos
        for s in range(Nscales):
            erode_its = np.ceil(1.0 / RESOLUTIONS[s]).astype("int")
            for z in range(Nslices):
                # M_ERODED = scipy.ndimage.binary_erosion(Ms[s][:, :, z] > .5, iterations=erode_its)
                for c in range(3):
                    Is[s][:, :, z, c] = (
                        grad2d(Is[s][:, :, z, c]) / 255.0
                    )  # * M_ERODED

    ########################################################

    print("Center the centers of gravity in the origin")

    if ref_type == "surface":
        cog_mesh_ras = np.mean(Pmesh, axis=0)
        Pmesh -= cog_mesh_ras
        REFaff[:-1, -1] = REFaff[:-1, -1] - cog_mesh_ras

    else:
        idx = np.where(REF > 0.1)
        cog_mri_vox = np.array(
            [[np.mean(idx[0])], [np.mean(idx[1])], [np.mean(idx[2])]]
        )
        cog_mri_ras = vox2ras(cog_mri_vox, REFaff)
        REFaff[:-1, -1] = REFaff[:-1, -1] - np.squeeze(cog_mri_ras)

    idx = np.where(Ms[-1] > 0)
    cog_photo_vox = np.array(
        [[np.mean(idx[0])], [np.mean(idx[1])], [np.mean(idx[2])]]
    )
    cog_photo_ras = vox2ras(cog_photo_vox, Affs[-1])
    for s in np.arange(Nscales):
        Affs[s][:-1, -1] = Affs[s][:-1, -1] - np.squeeze(cog_photo_ras)

    ########################################################

    print("Optimization")

    # Initialize
    t = None
    theta = None
    shear = None
    scaling = None
    sz = None
    t_reference = None
    theta_reference = None
    s_reference = None

    # Go over resolutions / modes
    n_modes = 2
    print("We will be running 2 modes: rigid, and affine")

    if options.rigid_only_for_photos:
        n_modes = 1
        print("We will only be running 1 mode only: rigid for everything")


    for mode_idx in range(n_modes):

        ref_type_iteration = ref_type
        k_dice_mri_iteration = K_DICE_MRI
        k_ncc_intermodality_iteration = K_NCC_INTERMODALITY
        k_surface_term_iteration = K_SURFACE_TERM
        allow_nonlin = False

        if mode_idx == 0:

            print("########################################################")
            print("###    First pass: no scaling / shearing allowed     ###")
            print("########################################################")

            allow_scaling_and_shear = False

        else:

            print("########################################################")
            print("###    Second pass: scaling / shearing is allowed    ###")
            print("########################################################")

            allow_scaling_and_shear = True


        for res in range(len(RESOLUTIONS)):

            print(
                "Working on resolution %d of %d (%.2f mm): %d iterations "
                % (res + 1, len(RESOLUTIONS), RESOLUTIONS[res], STEPS[res])
            )

            if ref_type == "surface":
                allow_s_reference = False
                ref_surface = Pmesh
            else:
                ref_surface = None

            volres = np.sqrt(np.sum(REFaff[:, :-1] ** 2, axis=0))
            sigmas = 0.5 * RESOLUTIONS[res] / volres
            REFsmooth = scipy.ndimage.gaussian_filter(REF, sigmas)
            allow_s_reference = (ref_type == "soft_mask")

            model = PhotoAligner(
                Is[res],
                Ms[res],
                Affs[res],
                REFsmooth,
                REFaff,
                ref_surface,
                t_ini=t,
                theta_ini=theta,
                shear_ini=shear,
                scaling_ini=scaling,
                sz_ini=sz,
                allow_sz=allow_z_stretch,
                t_reference_ini=t_reference,
                theta_reference_ini=theta_reference,
                s_reference_ini=s_reference,
                allow_s_reference=allow_s_reference,
                ref_type=ref_type_iteration,
                allow_nonlin=allow_nonlin,
                k_dice_mri=k_dice_mri_iteration,
                k_ncc_intermodality=k_ncc_intermodality_iteration,
                k_surface_term=k_surface_term_iteration,
                k_dice_slices=K_DICE_SLICES,
                k_ncc_slices=K_NCC_SLICES,
                k_regularizer=K_REGULARIZER,
                pad_ignore=PAD_AP,
                device=device,
                allow_scaling_and_shear=allow_scaling_and_shear
            )

            if alternate_optimization:
                optimizer2d = FullBatchLBFGS(model.parameters2d())
                optimizer3d = FullBatchLBFGS(model.parameters3d())
            else:
                optimizer = FullBatchLBFGS(model.parameters())

            loss_old = 1e10

            trigger_times = 0
            for epoch in range(STEPS[res]):

                # Compute loss with forward pass
                loss = model()[0]

                # optimize with BFGS
                if alternate_optimization:

                    def closure2d():
                        optimizer2d.zero_grad()
                        loss = model()[0]
                        return loss

                    def closure3d():
                        optimizer3d.zero_grad()
                        loss = model()[0]
                        return loss

                else:

                    def closure():
                        optimizer.zero_grad()
                        loss = model()[0]
                        return loss

                # optimizer.step(closure)
                if epoch == 1:
                    loss.backward()

                if alternate_optimization:
                    options2d = {
                        "closure": closure2d,
                        "current_loss": loss,
                        "max_ls": 75,
                    }
                    options3d = {
                        "closure": closure3d,
                        "current_loss": loss,
                        "max_ls": 75,
                    }

                    if epoch % 10 < 5:
                        (
                            loss,
                            _,
                            lr,
                            _,
                            F_eval,
                            G_eval,
                            _,
                            fail_flag,
                        ) = optimizer2d.step(options2d)
                    else:
                        (
                            loss,
                            _,
                            lr,
                            _,
                            F_eval,
                            G_eval,
                            _,
                            fail_flag,
                        ) = optimizer3d.step(options3d)
                else:
                    options = {
                        "closure": closure,
                        "current_loss": loss,
                        "max_ls": 75,
                    }
                    loss, _, lr, _, F_eval, G_eval, _, fail_flag = optimizer.step(
                        options
                    )

                if fail_flag:
                    print("Line search failed")
                    break

                # print step info
                loss = loss.cpu().detach().numpy()
                print("   Step %d, loss = %.10f" % (epoch + 1, loss), flush=True)

                if ((loss_old - loss) < TOL):
                    trigger_times += 1

                    if trigger_times >= 25:
                        print(
                            "   Decrease in loss below tolerance limit for the last 25 steps"
                        )
                        break
                else:
                    trigger_times = 0

                loss_old = loss

            # Retrieve model parameters
            t = model.t.cpu().detach().numpy()
            theta = model.theta.cpu().detach().numpy()
            shear = model.shear.cpu().detach().numpy()
            scaling = model.scaling.cpu().detach().numpy()
            sz = model.sz.cpu().detach().numpy()
            t_reference = model.t_reference.cpu().detach().numpy()
            theta_reference = model.theta_reference.cpu().detach().numpy()
            s_reference = model.s_reference.cpu().detach().numpy()


            # In the last resolution level, retrieve results before deleting the model
            if res == (len(RESOLUTIONS) - 1) and mode_idx == n_modes - 1:

                model.photo_vol = torch.Tensor(Is_copy).to(
                    device
                )  # TODO: I had a division by np.max(Is_copy) but got rid of it...
                model.photo_rearranged = torch.unsqueeze(
                    model.photo_vol.permute(3, 0, 1, 2), dim=0
                ).to(model.device)

                if ref_type == "surface":
                    (
                        _,
                        photo_resampled,
                        photo_aff,
                        mri_aff_combined,
                        Rt,
                        TvoxPhotos,
                    ) = model()
                    Rt = Rt.cpu().detach().numpy()
                else:
                    (
                        _,
                        photo_resampled,
                        photo_aff,
                        mri_aff_combined,
                        _,
                        TvoxPhotos,
                    ) = model()

                TvoxPhotos = TvoxPhotos.cpu().detach().numpy()
                mri_aff_combined = mri_aff_combined.cpu().detach().numpy()
                photo_resampled = photo_resampled.cpu().detach().numpy()
                photo_aff = photo_aff.cpu().detach().numpy()

            # Free up memory
            if alternate_optimization:
                del optimizer2d
                del optimizer3d
            else:
                del optimizer
            del model

    ########################################################

    print("Writing results to disk")

    if ref_type == "surface":

        MRIwrite(photo_resampled, photo_aff, output_photo_recon)

        Pmesh_rotated = np.matmul(
            np.concatenate([Pmesh, np.ones([Pmesh.shape[0], 1])], axis=1),
            Rt.transpose(),
        )[:, :-1]
        nib.freesurfer.write_geometry(
            output_registered_reference, Pmesh_rotated, Tmesh, volume_info=meta_mesh
        )

        reg_mask = output_directory + "registered_reference.mgz"
        MRIwrite(REF_orig, mri_aff_combined, reg_mask)

        print(
            "freeview -v %s -v %s -f %s"
            % (
                output_photo_recon,
                reg_mask,
                output_registered_reference,
            )
        )

    else:
        # Unless reference is soft, go back to original RAS space of reference before writing photo volume
        if ref_type == "soft_mask":
            MRIwrite(photo_resampled, photo_aff, output_photo_recon)
            reg_mask = output_directory + "registered_reference.mgz"
            MRIwrite(REF_orig, mri_aff_combined, reg_mask)
            print("freeview %s %s" % (output_photo_recon, reg_mask))

        else:
            T = np.matmul(mri_aff_combined, np.linalg.inv(REFaff_orig))
            Tinv = np.linalg.inv(T)
            MRIwrite(
                photo_resampled, np.matmul(Tinv, photo_aff), output_photo_recon
            )
            print("freeview %s %s" % (output_photo_recon, input_reference))


    print("All done!")


# ================================================================================================
#                                         Auxiliary functions and classes
# ================================================================================================
def seed_all(seed):
    # https://discuss.pytorch.org/t/reproducibility-with-all-the-bells-and-whistles/81097
    seed = 0 if not seed else seed
    print("[ Using Seed : ", seed, " ]")

    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.cuda.manual_seed(seed)
    np.random.seed(seed)
    np.random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

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

def grad2d(X, provide_gradients=False):
    h = np.array([-1, 0, 1])
    Gx = scipy.ndimage.convolve(X, np.reshape(h, [3, 1]))
    Gy = scipy.ndimage.convolve(X, np.reshape(h, [1, 3]))
    Gmodule = np.sqrt(Gx * Gx + Gy * Gy)
    if provide_gradients:
        return Gmodule, Gx, Gy
    else:
        return Gmodule

def MRIread(filename, dtype=None, im_only=False):

    assert filename.endswith(
        ('.nii', '.nii.gz', '.mgz')), 'Unknown data file: %s' % filename

    x = nib.load(filename)
    volume = x.get_data()
    aff = x.affine

    if dtype is not None:
        volume = volume.astype(dtype=dtype)

    if im_only:
        return volume
    else:
        return volume, aff

def MRIwrite(volume, aff, filename, dtype=None):

    if dtype is not None:
        volume = volume.astype(dtype=dtype)

    if aff is None:
        aff = np.eye(4)
    header = nib.Nifti1Header()
    nifty = nib.Nifti1Image(volume, aff, header)

    nib.save(nifty, filename)

def vox2ras(vox, vox2ras):
    vox2 = np.concatenate([vox, np.ones(shape=[1, vox.shape[1]])], axis=0)
    ras = np.matmul(vox2ras, vox2)[:-1, :]
    return ras

def grad3d(X, provide_gradients=False):
    h = np.array([-1, 0, 1])
    Gx = scipy.ndimage.convolve(X, np.reshape(h, [3, 1, 1]))
    Gy = scipy.ndimage.convolve(X, np.reshape(h, [1, 3, 1]))
    Gz = scipy.ndimage.convolve(X, np.reshape(h, [1, 1, 3]))
    Gmodule = np.sqrt(Gx * Gx + Gy * Gy + Gz * Gz)

    if provide_gradients:
        return Gmodule, Gx, Gy, Gz
    else:
        return Gmodule


class PhotoAligner(nn.Module):
    """
    Main class to perform alignment
    """

    def __init__(
        self,
        photo_vol,
        photo_mask_vol,
        photo_aff,
        mri_vol,
        mri_aff,
        ref_surf,
        device="cpu",
        allow_scaling_and_shear=False,
        allow_sz=False,
        k_dice_mri=0.8,
        k_ncc_intermodality=0.8,
        k_surface_term=0.8,
        k_dice_slices=0.1,
        k_ncc_slices=0.1,
        k_regularizer=0.01,
        pad_ignore=None,
        t_ini=None,
        theta_ini=None,
        shear_ini=None,
        scaling_ini=None,
        sz_ini=None,
        ref_type="mask",
        allow_s_reference=False,
        t_reference_ini=None,
        theta_reference_ini=None,
        s_reference_ini=None,
        allow_nonlin=False,
        field_ini=None,
        k_nonlinear=0.1,
    ):

        super().__init__()

        self.device = device
        self.ref_type = ref_type
        self.photo_vol = torch.Tensor(photo_vol).to(self.device)
        self.photo_rearranged = torch.unsqueeze(
            self.photo_vol.permute(3, 0, 1, 2), dim=0
        ).to(self.device)
        self.photo_mask_vol = torch.Tensor(photo_mask_vol).to(self.device)
        self.mask_rearranged = torch.unsqueeze(
            torch.unsqueeze(self.photo_mask_vol, dim=0), dim=0
        ).to(self.device)
        self.photo_aff = torch.Tensor(photo_aff).to(self.device)
        if ref_type == "surface":
            self.ref_surf = torch.Tensor(ref_surf).to(self.device)
        else:
            self.ref_surf = None
        self.mri_vol = torch.Tensor(mri_vol).to(self.device)
        self.mri_rearranged = torch.unsqueeze(
            torch.unsqueeze(self.mri_vol, dim=0), dim=0
        ).to(self.device)
        self.mri_aff = torch.Tensor(mri_aff).to(self.device)
        self.allow_scaling_and_shear = allow_scaling_and_shear
        self.photo_siz = photo_mask_vol.shape[:-1]
        self.Nslices = photo_mask_vol.shape[-1]
        self.k_dice_mri = k_dice_mri
        self.k_dice_slices = k_dice_slices
        self.k_ncc_slices = k_ncc_slices
        self.k_regularizer = k_regularizer
        self.k_ncc_intermodality = k_ncc_intermodality
        self.k_surface_term = k_surface_term
        self.allow_nonlin = allow_nonlin
        self.k_nonlinear = k_nonlinear

        if pad_ignore is not None:
            self.pad_ignore = [pad_ignore, pad_ignore]
        else:
            # Discover automatically
            idx = np.argwhere(
                np.sum(np.sum(photo_mask_vol, axis=0), axis=0) > 0
            )
            self.pad_ignore = [
                np.min(idx),
                photo_mask_vol.shape[-1] - 1 - np.max(idx),
            ]

        if t_ini is not None:
            self.t = torch.nn.Parameter(torch.tensor(t_ini).to(self.device))
        else:
            self.t = torch.nn.Parameter(
                torch.zeros(2, self.Nslices).to(self.device)
            )
        self.t.requires_grad = True

        if theta_ini is not None:
            self.theta = torch.nn.Parameter(
                torch.tensor(theta_ini).to(self.device)
            )
        else:
            self.theta = torch.nn.Parameter(
                torch.zeros(self.Nslices).to(self.device)
            )
        self.theta.requires_grad = True

        if allow_scaling_and_shear:
            if shear_ini is not None:
                self.shear = torch.nn.Parameter(
                    torch.tensor(shear_ini).to(self.device)
                )
            else:
                self.shear = torch.nn.Parameter(
                    torch.zeros(2, self.Nslices).to(self.device)
                )
            self.shear.requires_grad = True
            if scaling_ini is not None:
                self.scaling = torch.nn.Parameter(
                    torch.tensor(scaling_ini).to(self.device)
                )
            else:
                self.scaling = torch.nn.Parameter(
                    torch.zeros(2, self.Nslices).to(self.device)
                )
            self.scaling.requires_grad = True
        else:
            if shear_ini is not None:
                self.shear = torch.tensor(shear_ini).to(self.device)
            else:
                self.shear = torch.zeros(2, self.Nslices).to(self.device)
            if scaling_ini is not None:
                self.scaling = torch.tensor(scaling_ini).to(self.device)
            else:
                self.scaling = torch.zeros(2, self.Nslices).to(self.device)

        if allow_sz:
            if sz_ini is not None:
                self.sz = torch.nn.Parameter(
                    torch.tensor(sz_ini).to(self.device)
                )
            else:
                self.sz = torch.nn.Parameter(torch.zeros(1).to(self.device))
            self.sz.requires_grad = True
        else:
            if sz_ini is not None:
                self.sz = torch.tensor(sz_ini).to(self.device)
            else:
                self.sz = torch.zeros(1).to(self.device)

        if allow_s_reference:
            if s_reference_ini is not None:
                self.s_reference = torch.nn.Parameter(
                    torch.tensor(s_reference_ini).to(self.device)
                )
            else:
                self.s_reference = torch.nn.Parameter(
                    torch.zeros(1).to(self.device)
                )
            self.s_reference.requires_grad = True
        else:
            if s_reference_ini is not None:
                self.s_reference = torch.tensor(s_reference_ini).to(self.device)
            else:
                self.s_reference = torch.zeros(1).to(self.device)

        if t_reference_ini is not None:
            self.t_reference = torch.nn.Parameter(
                torch.tensor(t_reference_ini).to(self.device)
            )
        else:
            self.t_reference = torch.nn.Parameter(
                torch.zeros(3).to(self.device)
            )
        self.t_reference.requires_grad = True

        if theta_reference_ini is not None:
            self.theta_reference = torch.nn.Parameter(
                torch.tensor(theta_reference_ini).to(self.device)
            )
        else:
            self.theta_reference = torch.nn.Parameter(
                torch.zeros(3).to(self.device)
            )
        self.theta_reference.requires_grad = True

        if allow_nonlin:
            if field_ini is not None:
                self.field = torch.nn.Parameter(
                    torch.tensor(field_ini).to(self.device)
                )
            else:
                self.field = torch.nn.Parameter(
                    torch.zeros(5, 5, 2, self.Nslices).to(self.device)
                )
            self.field.requires_grad = True
        else:
            if field_ini is not None:
                self.field = torch.tensor(field_ini).to(self.device)
            else:
                self.field = None

        # create sampling grid
        vectors = [torch.arange(0, s) for s in self.photo_mask_vol.shape]
        self.grids = torch.stack(torch.meshgrid(vectors)).to(self.device)

        self.DELTA = 0.1 / np.log(photo_vol.shape[1])
        self.pad_ap = 3

    def forward(self):

        # We scale angles / shearings / scalings as a simple form of preconditioning (which shouldn't be needed with bfgs, but whatever...)

        # Parameters of photos
        theta_f = (
            self.theta / 180 * torch.tensor(np.pi)
        )  # + 0.1  # degrees -> radians
        shear_f = self.shear / 100  # + 0.1 # percentages
        scaling_f = torch.exp(
            self.scaling / 20
        )  # ensures positive and symmetry around 1 in log scale
        t_f = self.t + self.DELTA  # no scaling
        sz_f = torch.exp(
            self.sz / 20
        )  # ensures positive and symmetry around 1 in log scale

        # Parameters of reference volume, if it exists
        theta_reference_f = (
            self.theta_reference / 180 * torch.tensor(np.pi)
        )  # degrees -> radians
        t_reference_f = self.t_reference  # no scaling
        s_reference_f = torch.exp(
            self.s_reference / 20
        )  # ensures positive and symmetry around 1 in log scale
        if self.field is not None:
            # make it a % of the dimensions
            field_x = self.field[0, :, :, :] * torch.Tensor(
                np.array(self.photo_siz[0] / 100.0)
            )
            field_y = self.field[1, :, :, :] * torch.Tensor(
                np.array(self.photo_siz[1] / 100.0)
            )
            field_pixels = torch.stack([field_x, field_y])
        else:
            field_pixels = None

        if (
            torch.any(torch.isnan(theta_f))
            or torch.any(torch.isnan(shear_f))
            or torch.any(torch.isnan(scaling_f))
            or torch.any(torch.isnan(t_f))
            or torch.any(torch.isnan(sz_f))
            or torch.any(torch.isnan(theta_reference_f))
            or torch.any(torch.isnan(t_reference_f))
            or torch.any(torch.isnan(s_reference_f))
        ):
            kk = 1

        # Prepare 2D matrices for the photos
        M = torch.zeros([4, 4, self.Nslices]).to(self.device)
        M[0, 0, :] = scaling_f[0, :] * (
            torch.cos(theta_f) - shear_f[0, :] * torch.sin(theta_f)
        )
        M[0, 2, :] = scaling_f[0, :] * (
            shear_f[1, :] * torch.cos(theta_f)
            - (1 + shear_f[0, :] * shear_f[1, :]) * torch.sin(theta_f)
        )
        M[0, 3, :] = t_f[0, :]
        M[1, 1, :] = 1
        M[2, 0, :] = scaling_f[1, :] * (
            torch.sin(theta_f) + shear_f[0, :] * torch.cos(theta_f)
        )
        M[2, 2, :] = scaling_f[1, :] * (
            shear_f[1, :] * torch.sin(theta_f)
            + (1 + shear_f[0, :] * shear_f[1, :]) * torch.cos(theta_f)
        )
        M[2, 3, :] = t_f[1, :]
        M[3, 3, :] = 1

        # update mesh grids for photos
        # First, upscale field
        if self.field is not None:

            field_fullsiz = torch.nn.Upsample(
                size=self.photo_vol.shape[:-2],
                align_corners=True,
                mode="bilinear",
            )(field_pixels)
            field_fullsiz_rearranged = field_fullsiz.permute(0, 2, 3, 1)

            # cost_field = torch.mean(torch.sqrt(torch.clamp(field_pixels * field_pixels, min=1e-5)))
            grad_xx = (
                field_fullsiz[0, 2:, 1:-1, :] - field_fullsiz[0, :-2, 1:-1, :]
            ) / 2.0
            grad_xy = (
                field_fullsiz[0, 1:-1, 2:, :] - field_fullsiz[0, 1:-1, :-2, :]
            ) / 2.0
            grad_x = torch.sqrt(
                torch.clamp(grad_xx * grad_xx + grad_xy * grad_xy, min=1e-5)
            )
            grad_yx = (
                field_fullsiz[1, 2:, 1:-1, :] - field_fullsiz[1, :-2, 1:-1, :]
            ) / 2.0
            grad_yy = (
                field_fullsiz[1, 1:-1, 2:, :] - field_fullsiz[1, 1:-1, :-2, :]
            ) / 2.0
            grad_y = torch.sqrt(
                torch.clamp(grad_yx * grad_yx + grad_yy * grad_yy, min=1e-5)
            )
            cost_field = torch.mean(grad_x) + torch.mean(grad_y)

        else:
            field_fullsiz = None
            field_fullsiz_rearranged = None
            cost_field = torch.zeros(1).to(self.device)

        # update mesh grids for photos
        grids_new = torch.zeros(self.grids.shape).to(self.device)
        photo_aff = torch.zeros(4, 4).to(self.device)
        photo_aff[:, :] = self.photo_aff
        photo_aff[1, 2] = self.photo_aff[1, 2] * sz_f
        T = torch.zeros([4, 4, self.Nslices]).to(self.device)
        for z in range(self.Nslices):
            T[:, :, z] = torch.matmul(
                torch.matmul(torch.inverse(photo_aff), M[:, :, z]), photo_aff
            )
            if field_fullsiz_rearranged is None:
                for d in range(3):
                    grids_new[d, :, :, z] = (
                        T[d, 0, z] * self.grids[0, :, :, z]
                        + T[d, 1, z] * self.grids[1, :, :, z]
                        + T[d, 2, z] * self.grids[2, :, :, z]
                        + T[d, 3, z]
                    )
            else:
                for d in range(2):
                    grids_new[d, :, :, z] = (
                        T[d, 0, z] * self.grids[0, :, :, z]
                        + T[d, 1, z] * self.grids[1, :, :, z]
                        + T[d, 2, z] * self.grids[2, :, :, z]
                        + T[d, 3, z]
                        + field_fullsiz_rearranged[d, :, :, z]
                    )
                grids_new[2, :, :, z] = (
                    T[2, 0, z] * self.grids[0, :, :, z]
                    + T[2, 1, z] * self.grids[1, :, :, z]
                    + T[2, 2, z] * self.grids[2, :, :, z]
                    + T[2, 3, z]
                )

        # Resample photos and masks
        # We need to make the new grid compatible with grid_resample...
        grids_new = torch.unsqueeze(grids_new, 0)
        grids_new = grids_new.permute(0, 2, 3, 4, 1)
        for i in range(3):
            grids_new[:, :, :, :, i] = 2 * (
                grids_new[:, :, :, :, i] / (self.photo_mask_vol.shape[i] - 1)
                - 0.5
            )
        # Not sure why, but channels need to be reversed
        grids_new = grids_new[..., [2, 1, 0]]
        photo_resampled = nnf.grid_sample(
            self.photo_rearranged,
            grids_new,
            align_corners=True,
            mode="bilinear",
            padding_mode="zeros",
        )
        photo_resampled = torch.squeeze(photo_resampled.permute(2, 3, 4, 1, 0))
        mask_resampled = nnf.grid_sample(
            self.mask_rearranged,
            grids_new,
            align_corners=True,
            mode="bilinear",
            padding_mode="zeros",
        )
        mask_resampled = torch.squeeze(mask_resampled)

        # Now work on the reference
        Rx = torch.zeros([4, 4]).to(self.device)
        Rx[0, 0] = 1
        Rx[1, 1] = torch.cos(theta_reference_f[0])
        Rx[1, 2] = -torch.sin(theta_reference_f[0])
        Rx[2, 1] = torch.sin(theta_reference_f[0])
        Rx[2, 2] = torch.cos(theta_reference_f[0])
        Rx[3, 3] = 1

        Ry = torch.zeros([4, 4]).to(self.device)
        Ry[0, 0] = torch.cos(theta_reference_f[1])
        Ry[0, 2] = torch.sin(theta_reference_f[1])
        Ry[1, 1] = 1
        Ry[2, 0] = -torch.sin(theta_reference_f[1])
        Ry[2, 2] = torch.cos(theta_reference_f[1])
        Ry[3, 3] = 1

        Rz = torch.zeros([4, 4]).to(self.device)
        Rz[0, 0] = torch.cos(theta_reference_f[2])
        Rz[0, 1] = -torch.sin(theta_reference_f[2])
        Rz[1, 0] = torch.sin(theta_reference_f[2])
        Rz[1, 1] = torch.cos(theta_reference_f[2])
        Rz[2, 2] = 1
        Rz[3, 3] = 1

        trans_and_scale = torch.eye(4).to(self.device)
        trans_and_scale[:-1, -1] = t_reference_f
        trans_and_scale[0, 0] = s_reference_f
        trans_and_scale[1, 1] = s_reference_f
        trans_and_scale[2, 2] = s_reference_f

        Rt = torch.matmul(
            trans_and_scale, torch.matmul(torch.matmul(Rx, Ry), Rz)
        )

        if self.ref_type == "surface":

            # Combine transformation of surface with ras2vox
            combo = torch.matmul(torch.inverse(photo_aff), Rt)
            xx = (
                combo[0, 0] * self.ref_surf[:, 0]
                + combo[0, 1] * self.ref_surf[:, 1]
                + combo[0, 2] * self.ref_surf[:, 2]
                + combo[0, 3]
            )
            yy = (
                combo[1, 0] * self.ref_surf[:, 0]
                + combo[1, 1] * self.ref_surf[:, 1]
                + combo[1, 2] * self.ref_surf[:, 2]
                + combo[1, 3]
            )
            zz = (
                combo[2, 0] * self.ref_surf[:, 0]
                + combo[2, 1] * self.ref_surf[:, 1]
                + combo[2, 2] * self.ref_surf[:, 2]
                + combo[2, 3]
            )

            ok = (
                (xx > 0)
                & (yy > 0)
                & (zz > 0)
                & (xx < photo_resampled.shape[0] - 1)
                & (yy < photo_resampled.shape[1] - 1)
                & (zz < photo_resampled.shape[2] - 1)
            )

            x = xx[ok]
            y = yy[ok]
            z = zz[ok]

            # Interpolate in a nice, differentiable way. We follow https://en.wikipedia.org/wiki/Trilinear_interpolation
            xf = torch.floor(x).long()
            yf = torch.floor(y).long()
            zf = torch.floor(z).long()
            xd = (x - xf).unsqueeze(1)
            yd = (y - yf).unsqueeze(1)
            zd = (z - zf).unsqueeze(1)

            c000 = photo_resampled[xf, yf, zf, :]
            c001 = photo_resampled[xf, yf, zf + 1, :]
            c010 = photo_resampled[xf, yf + 1, zf, :]
            c011 = photo_resampled[xf, yf + 1, zf + 1, :]
            c100 = photo_resampled[xf + 1, yf, zf, :]
            c101 = photo_resampled[xf + 1, yf, zf + 1, :]
            c110 = photo_resampled[xf + 1, yf + 1, zf, :]
            c111 = photo_resampled[xf + 1, yf + 1, zf + 1, :]

            c00 = c000 * (1 - xd) + c100 * xd
            c01 = c001 * (1 - xd) + c101 * xd
            c10 = c010 * (1 - xd) + c110 * xd
            c11 = c011 * (1 - xd) + c111 * xd

            c0 = c00 * (1 - yd) + c10 * yd
            c1 = c01 * (1 - yd) + c11 * yd

            c = c0 * (1 - zd) + c1 * zd

            av_grad = torch.sum(c) / torch.numel(ok) / 3.0

        # Now let's work on the deformed volume
        grids_new_mri = torch.zeros(self.grids.shape).to(self.device)
        mri_aff_combined = torch.matmul(Rt, self.mri_aff)
        D = torch.matmul(torch.inverse(mri_aff_combined), photo_aff)

        for d in range(3):
            grids_new_mri[d, :, :, :] = (
                D[d, 0] * self.grids[0, :, :, :]
                + D[d, 1] * self.grids[1, :, :, :]
                + D[d, 2] * self.grids[2, :, :, :]
                + grids_new_mri[d, :, :, :]
                + D[d, 3]
            )

        grids_new_mri = torch.unsqueeze(grids_new_mri, 0)
        grids_new_mri = grids_new_mri.permute(0, 2, 3, 4, 1)
        for i in range(3):
            grids_new_mri[:, :, :, :, i] = 2 * (
                grids_new_mri[:, :, :, :, i] / (self.mri_vol.shape[i] - 1) - 0.5
            )
        # Not sure why, but channels need to be reversed
        grids_new_mri = grids_new_mri[..., [2, 1, 0]]
        mri_resampled = nnf.grid_sample(
            self.mri_rearranged,
            grids_new_mri,
            align_corners=True,
            mode="bilinear",
            padding_mode="zeros",
        )
        mri_resampled = torch.squeeze(mri_resampled)

        if self.ref_type == "image":

            # I now do slice-wise ncc, rather than slice-wise lncc
            if False:

                def slice_wise_LNCC(I, J, win_semilen=4):

                    # set window size
                    win = [1 + 2 * win_semilen, 1 + 2 * win_semilen, 1]
                    sum_filt = torch.ones([1, 1, *win]).to(I.device)
                    pad_no = win_semilen
                    stride = (1, 1, 1)
                    padding = (pad_no, pad_no, 0)

                    # get convolution function
                    conv = getattr(torch.nn.functional, "conv3d")

                    Ipad = I[None, None, ...] / torch.max(I)
                    Jpad = J[None, None, ...] / torch.max(J)

                    # compute CC squares
                    I2 = Ipad * Ipad
                    J2 = Jpad * Jpad
                    IJ = Ipad * Jpad

                    I_sum = conv(Ipad, sum_filt, stride=stride, padding=padding)
                    J_sum = conv(Jpad, sum_filt, stride=stride, padding=padding)
                    I2_sum = conv(I2, sum_filt, stride=stride, padding=padding)
                    J2_sum = conv(J2, sum_filt, stride=stride, padding=padding)
                    IJ_sum = conv(IJ, sum_filt, stride=stride, padding=padding)

                    win_size = np.prod(win)
                    u_I = I_sum / win_size
                    u_J = J_sum / win_size

                    cross = (
                        IJ_sum
                        - u_J * I_sum
                        - u_I * J_sum
                        + u_I * u_J * win_size
                    )
                    I_var = I2_sum - 2 * u_I * I_sum + u_I * u_I * win_size
                    J_var = J2_sum - 2 * u_J * J_sum + u_J * u_J * win_size

                    cc = (cross * cross) / (I_var * J_var + 1e-6)

                    return torch.mean(cc)

                ncc_mri = (
                    slice_wise_LNCC(
                        mri_resampled[
                            :, :, self.pad_ignore[0] : -self.pad_ignore[1]
                        ],
                        photo_resampled[
                            :, :, self.pad_ignore[0] : -self.pad_ignore[1], 0
                        ],
                    )
                    / 3.0
                    + slice_wise_LNCC(
                        mri_resampled[
                            :, :, self.pad_ignore[0] : -self.pad_ignore[1]
                        ],
                        photo_resampled[
                            :, :, self.pad_ignore[0] : -self.pad_ignore[1], 1
                        ],
                    )
                    / 3.0
                    + slice_wise_LNCC(
                        mri_resampled[
                            :, :, self.pad_ignore[0] : -self.pad_ignore[1]
                        ],
                        photo_resampled[
                            :, :, self.pad_ignore[0] : -self.pad_ignore[1], 2
                        ],
                    )
                    / 3.0
                )

            else:
                nccs = torch.zeros(
                    3, self.Nslices - self.pad_ignore[0] - self.pad_ignore[1]
                ).to(self.device)
                for z in range(
                    self.Nslices - self.pad_ignore[0] - self.pad_ignore[1]
                ):
                    x = mri_resampled[:, :, z + self.pad_ignore[0]]
                    mx = torch.mean(x)
                    vx = x - mx
                    for c in range(3):
                        y = photo_resampled[:, :, z + self.pad_ignore[0], c]
                        my = torch.mean(y)
                        vy = y - my
                        nccs[c, z] = (
                            torch.mean(vx * vy)
                            / torch.sqrt(
                                torch.clamp(torch.mean(vx**2), min=1e-5)
                            )
                            / torch.sqrt(
                                torch.clamp(torch.mean(vy**2), min=1e-5)
                            )
                        )

                ncc_mri = torch.mean(nccs)

        else:
            num = torch.sum(2 * (mri_resampled * mask_resampled))
            den = torch.clamp(
                torch.sum(mri_resampled + mask_resampled), min=1e-5
            )
            dice_mri = num / den

        dices_slices = torch.zeros(
            self.Nslices - 1 - self.pad_ignore[0] - self.pad_ignore[1]
        ).to(self.device)
        for z in range(
            self.Nslices - 1 - self.pad_ignore[0] - self.pad_ignore[1]
        ):
            dices_slices[z] = torch.sum(
                2
                * (
                    mask_resampled[:, :, z + self.pad_ignore[0]]
                    * mask_resampled[:, :, z + self.pad_ignore[0] + 1]
                )
            ) / torch.clamp(
                torch.sum(
                    mask_resampled[:, :, z + self.pad_ignore[0]]
                    + mask_resampled[:, :, z + self.pad_ignore[0] + 1]
                ),
                min=1e-5,
            )
        dice_slices = torch.mean(dices_slices)

        nccs_slices = torch.zeros(
            self.Nslices - 1 - self.pad_ignore[0] - self.pad_ignore[1]
        ).to(self.device)
        for z in range(
            self.Nslices - 1 - self.pad_ignore[0] - self.pad_ignore[1]
        ):
            x = photo_resampled[:, :, z + self.pad_ignore[0], ...]
            y = photo_resampled[:, :, z + self.pad_ignore[0] + 1, ...]
            m = (
                mask_resampled[:, :, z + self.pad_ignore[0]]
                * mask_resampled[:, :, z + self.pad_ignore[0] + 1]
            )
            if len(x.shape) == 2:  # grayscale
                m3 = m
            else:  # rgb
                m3 = m.repeat([3, 1, 1]).permute([1, 2, 0])
            bottom = torch.clamp(torch.sum(m3), min=1e-5)
            mx = torch.sum(x * m3) / bottom
            my = torch.sum(y * m3) / bottom
            vx = x - mx
            vy = y - my
            nccs_slices[z] = (
                torch.sum(vx * vy * m3)
                / torch.sqrt(torch.clamp(torch.sum(m3 * (vx**2)), min=1e-5))
                / torch.sqrt(torch.clamp(torch.sum(m3 * (vy**2)), min=1e-5))
            )
        ncc_slices = torch.mean(nccs_slices)

        aff_regularizers = torch.abs(torch.sum(self.scaling, dim=0) / 20)
        regularizer = torch.mean(aff_regularizers)

        loss_photos = (
            -self.k_dice_slices * dice_slices
            - self.k_ncc_slices * ncc_slices
            + self.k_regularizer * regularizer
            + self.k_nonlinear * cost_field
        )

        if self.ref_type == "mask" or self.ref_type == "soft_mask":
            loss = loss_photos - self.k_dice_mri * dice_mri

        elif self.ref_type == "image":
            loss = loss_photos - self.k_ncc_intermodality * ncc_mri

        else:
            loss = (
                loss_photos
                - self.k_surface_term * av_grad
                - self.k_dice_mri * dice_mri
            )
            # loss = loss - 0.01 * torch.abs(torch.log(sz_f))

        # loss = loss + 0.001 * (torch.sum(torch.square((self.t[:, 0:self.pad_ap] - self.DELTA))))
        # loss = loss + 0.001 * (torch.sum(torch.square((self.t[:, -self.pad_ap:] - self.DELTA))))

        # loss = loss + 0.001 * (torch.sum(torch.square((self.theta[0:self.pad_ap]))))
        # loss = loss + 0.001 * (torch.sum(torch.square((self.theta[-self.pad_ap:]))))

        # loss = loss + 0.001 * (torch.sum(torch.square((self.shear[:, 0:self.pad_ap]))))
        # loss = loss + 0.001 * (torch.sum(torch.square((self.shear[:, -self.pad_ap:]))))

        # loss = loss + 0.001 * (torch.sum(torch.square((self.scaling[:, 0:self.pad_ap]))))
        # loss = loss + 0.001 * (torch.sum(torch.square((self.scaling[:, -self.pad_ap:]))))

        TINY = 1e-6
        loss = loss + TINY * (torch.mean(torch.square((self.t - self.DELTA))))
        loss = loss + TINY * (torch.mean(torch.square((self.theta))))
        loss = loss + TINY * (torch.mean(torch.square((self.shear))))
        loss = loss + TINY * (torch.mean(torch.square((self.scaling))))

        if torch.isnan(loss):
            kk = 1

        return loss, photo_resampled, photo_aff, mri_aff_combined, Rt, T

    def parameters2d(self):
        yield self.t
        yield self.theta
        if isinstance(self.shear, nn.Parameter):
            yield self.shear
        if isinstance(self.scaling, nn.Parameter):
            yield self.scaling
        if isinstance(self.sz, nn.Parameter):
            yield self.sz
        if isinstance(self.field, nn.Parameter):
            yield self.field

    def parameters3d(self):
        yield self.t_reference
        yield self.theta_reference
        if isinstance(self.s_reference, nn.Parameter):
            yield self.s_reference


def is_legal(v):
    """
    Checks that tensor is not NaN or Inf.

    Inputs:
        v (tensor): tensor to be checked

    """
    legal = not torch.isnan(v).any() and not torch.isinf(v)

    return legal


def polyinterp(points, x_min_bound=None, x_max_bound=None, plot=False):
    """
    Gives the minimizer and minimum of the interpolating polynomial over given points
    based on function and derivative information. Defaults to bisection if no critical
    points are valid.

    Based on polyinterp.m Matlab function in minFunc by Mark Schmidt with some slight
    modifications.

    Implemented by: Hao-Jun Michael Shi and Dheevatsa Mudigere
    Last edited 12/6/18.

    Inputs:
        points (nparray): two-dimensional array with each point of form [x f g]
        x_min_bound (float): minimum value that brackets minimum (default: minimum of points)
        x_max_bound (float): maximum value that brackets minimum (default: maximum of points)
        plot (bool): plot interpolating polynomial

    Outputs:
        x_sol (float): minimizer of interpolating polynomial
        F_min (float): minimum of interpolating polynomial

    Note:
      . Set f or g to np.nan if they are unknown

    """
    no_points = points.shape[0]
    order = np.sum(1 - np.isnan(points[:, 1:3]).astype("int")) - 1

    x_min = np.min(points[:, 0])
    x_max = np.max(points[:, 0])

    # compute bounds of interpolation area
    if x_min_bound is None:
        x_min_bound = x_min
    if x_max_bound is None:
        x_max_bound = x_max

    # explicit formula for quadratic interpolation
    if no_points == 2 and order == 2 and plot is False:
        # Solution to quadratic interpolation is given by:
        # a = -(f1 - f2 - g1(x1 - x2))/(x1 - x2)^2
        # x_min = x1 - g1/(2a)
        # if x1 = 0, then is given by:
        # x_min = - (g1*x2^2)/(2(f2 - f1 - g1*x2))

        if points[0, 0] == 0:
            x_sol = (
                -points[0, 2]
                * points[1, 0] ** 2
                / (2 * (points[1, 1] - points[0, 1] - points[0, 2] * points[1, 0]))
            )
        else:
            a = (
                -(
                    points[0, 1]
                    - points[1, 1]
                    - points[0, 2] * (points[0, 0] - points[1, 0])
                )
                / (points[0, 0] - points[1, 0]) ** 2
            )
            x_sol = points[0, 0] - points[0, 2] / (2 * a)

        x_sol = np.minimum(np.maximum(x_min_bound, x_sol), x_max_bound)

    # explicit formula for cubic interpolation
    elif no_points == 2 and order == 3 and plot is False:
        # Solution to cubic interpolation is given by:
        # d1 = g1 + g2 - 3((f1 - f2)/(x1 - x2))
        # d2 = sqrt(d1^2 - g1*g2)
        # x_min = x2 - (x2 - x1)*((g2 + d2 - d1)/(g2 - g1 + 2*d2))
        d1 = (
            points[0, 2]
            + points[1, 2]
            - 3 * ((points[0, 1] - points[1, 1]) / (points[0, 0] - points[1, 0]))
        )
        d2 = np.sqrt(d1**2 - points[0, 2] * points[1, 2])
        if np.isreal(d2):
            x_sol = points[1, 0] - (points[1, 0] - points[0, 0]) * (
                (points[1, 2] + d2 - d1) / (points[1, 2] - points[0, 2] + 2 * d2)
            )
            x_sol = np.minimum(np.maximum(x_min_bound, x_sol), x_max_bound)
        else:
            x_sol = (x_max_bound + x_min_bound) / 2

    # solve linear system
    else:
        # define linear constraints
        A = np.zeros((0, order + 1))
        b = np.zeros((0, 1))

        # add linear constraints on function values
        for i in range(no_points):
            if not np.isnan(points[i, 1]):
                constraint = np.zeros((1, order + 1))
                for j in range(order, -1, -1):
                    constraint[0, order - j] = points[i, 0] ** j
                A = np.append(A, constraint, 0)
                b = np.append(b, points[i, 1])

        # add linear constraints on gradient values
        for i in range(no_points):
            if not np.isnan(points[i, 2]):
                constraint = np.zeros((1, order + 1))
                for j in range(order):
                    constraint[0, j] = (order - j) * points[i, 0] ** (order - j - 1)
                A = np.append(A, constraint, 0)
                b = np.append(b, points[i, 2])

        # check if system is solvable
        if A.shape[0] != A.shape[1] or np.linalg.matrix_rank(A) != A.shape[0]:
            x_sol = (x_min_bound + x_max_bound) / 2
            f_min = np.Inf
        else:
            # solve linear system for interpolating polynomial
            coeff = np.linalg.solve(A, b)

            # compute critical points
            dcoeff = np.zeros(order)
            for i in range(len(coeff) - 1):
                dcoeff[i] = coeff[i] * (order - i)

            crit_pts = np.array([x_min_bound, x_max_bound])
            crit_pts = np.append(crit_pts, points[:, 0])

            if not np.isinf(dcoeff).any():
                roots = np.roots(dcoeff)
                crit_pts = np.append(crit_pts, roots)

            # test critical points
            f_min = np.Inf
            x_sol = (x_min_bound + x_max_bound) / 2  # defaults to bisection
            for crit_pt in crit_pts:
                if (
                    np.isreal(crit_pt)
                    and crit_pt >= x_min_bound
                    and crit_pt <= x_max_bound
                ):
                    F_cp = np.polyval(coeff, crit_pt)
                    if np.isreal(F_cp) and F_cp < f_min:
                        x_sol = np.real(crit_pt)
                        f_min = np.real(F_cp)

            if plot:
                plt.figure()
                x = np.arange(
                    x_min_bound, x_max_bound, (x_max_bound - x_min_bound) / 10000
                )
                f = np.polyval(coeff, x)
                plt.plot(x, f)
                plt.plot(x_sol, f_min, "x")

    return x_sol


class LBFGS(Optimizer):
    """
    Implements the L-BFGS algorithm. Compatible with multi-batch and full-overlap
    L-BFGS implementations and (stochastic) Powell damping. Partly based on the
    original L-BFGS implementation in PyTorch, Mark Schmidt's minFunc MATLAB code,
    and Michael Overton's weak Wolfe line search MATLAB code.

    Implemented by: Hao-Jun Michael Shi and Dheevatsa Mudigere
    Last edited 10/20/20.

    Warnings:
      . Does not support per-parameter options and parameter groups.
      . All parameters have to be on a single device.

    Inputs:
        lr (float): steplength or learning rate (default: 1)
        history_size (int): update history size (default: 10)
        line_search (str): designates line search to use (default: 'Wolfe')
            Options:
                'None': uses steplength designated in algorithm
                'Armijo': uses Armijo backtracking line search
                'Wolfe': uses Armijo-Wolfe bracketing line search
        dtype: data type (default: torch.float)
        debug (bool): debugging mode

    References:
    [1] Berahas, Albert S., Jorge Nocedal, and Martin Takác. "A Multi-Batch L-BFGS
        Method for Machine Learning." Advances in Neural Information Processing
        Systems. 2016.
    [2] Bollapragada, Raghu, et al. "A Progressive Batching L-BFGS Method for Machine
        Learning." International Conference on Machine Learning. 2018.
    [3] Lewis, Adrian S., and Michael L. Overton. "Nonsmooth Optimization via Quasi-Newton
        Methods." Mathematical Programming 141.1-2 (2013): 135-163.
    [4] Liu, Dong C., and Jorge Nocedal. "On the Limited Memory BFGS Method for
        Large Scale Optimization." Mathematical Programming 45.1-3 (1989): 503-528.
    [5] Nocedal, Jorge. "Updating Quasi-Newton Matrices With Limited Storage."
        Mathematics of Computation 35.151 (1980): 773-782.
    [6] Nocedal, Jorge, and Stephen J. Wright. "Numerical Optimization." Springer New York,
        2006.
    [7] Schmidt, Mark. "minFunc: Unconstrained Differentiable Multivariate Optimization
        in Matlab." Software available at http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html
        (2005).
    [8] Schraudolph, Nicol N., Jin Yu, and Simon Günter. "A Stochastic Quasi-Newton
        Method for Online Convex Optimization." Artificial Intelligence and Statistics.
        2007.
    [9] Wang, Xiao, et al. "Stochastic Quasi-Newton Methods for Nonconvex Stochastic
        Optimization." SIAM Journal on Optimization 27.2 (2017): 927-956.

    """

    def __init__(
        self,
        params,
        lr=1.0,
        history_size=10,
        line_search="Wolfe",
        dtype=torch.float,
        debug=False,
    ):

        # ensure inputs are valid
        if not 0.0 <= lr:
            raise ValueError("Invalid learning rate: {}".format(lr))
        if not 0 <= history_size:
            raise ValueError("Invalid history size: {}".format(history_size))
        if line_search not in ["Armijo", "Wolfe", "None"]:
            raise ValueError("Invalid line search: {}".format(line_search))

        defaults = dict(
            lr=lr,
            history_size=history_size,
            line_search=line_search,
            dtype=dtype,
            debug=debug,
        )
        super(LBFGS, self).__init__(params, defaults)

        if len(self.param_groups) != 1:
            raise ValueError(
                "L-BFGS doesn't support per-parameter options " "(parameter groups)"
            )

        self._params = self.param_groups[0]["params"]
        self._numel_cache = None

        state = self.state["global_state"]
        state.setdefault("n_iter", 0)
        state.setdefault("curv_skips", 0)
        state.setdefault("fail_skips", 0)
        state.setdefault("H_diag", 1)
        state.setdefault("fail", True)

        state["old_dirs"] = []
        state["old_stps"] = []

    def _numel(self):
        if self._numel_cache is None:
            self._numel_cache = reduce(
                lambda total, p: total + p.numel(), self._params, 0
            )
        return self._numel_cache

    def _gather_flat_grad(self):
        views = []
        for p in self._params:
            if p.grad is None:
                view = p.data.new(p.data.numel()).zero_()
            elif p.grad.data.is_sparse:
                view = p.grad.data.to_dense().view(-1)
            else:
                view = p.grad.data.view(-1)
            views.append(view)
        return torch.cat(views, 0)

    def _add_update(self, step_size, update):
        offset = 0
        for p in self._params:
            numel = p.numel()
            # view as to avoid deprecated pointwise semantics
            p.data.add_(step_size, update[offset : offset + numel].view_as(p.data))
            offset += numel
        assert offset == self._numel()

    def _copy_params(self):
        current_params = []
        for param in self._params:
            current_params.append(deepcopy(param.data))
        return current_params

    def _load_params(self, current_params):
        i = 0
        for param in self._params:
            param.data[:] = current_params[i]
            i += 1

    def line_search(self, line_search):
        """
        Switches line search option.

        Inputs:
            line_search (str): designates line search to use
                Options:
                    'None': uses steplength designated in algorithm
                    'Armijo': uses Armijo backtracking line search
                    'Wolfe': uses Armijo-Wolfe bracketing line search

        """

        group = self.param_groups[0]
        group["line_search"] = line_search

        return

    def two_loop_recursion(self, vec):
        """
        Performs two-loop recursion on given vector to obtain Hv.

        Inputs:
            vec (tensor): 1-D tensor to apply two-loop recursion to

        Output:
            r (tensor): matrix-vector product Hv

        """

        group = self.param_groups[0]
        history_size = group["history_size"]

        state = self.state["global_state"]
        old_dirs = state.get("old_dirs")  # change in gradients
        old_stps = state.get("old_stps")  # change in iterates
        H_diag = state.get("H_diag")

        # compute the product of the inverse Hessian approximation and the gradient
        num_old = len(old_dirs)

        if "rho" not in state:
            state["rho"] = [None] * history_size
            state["alpha"] = [None] * history_size
        rho = state["rho"]
        alpha = state["alpha"]

        for i in range(num_old):
            rho[i] = 1.0 / old_stps[i].dot(old_dirs[i])

        q = vec
        for i in range(num_old - 1, -1, -1):
            alpha[i] = old_dirs[i].dot(q) * rho[i]
            q.add_(-alpha[i], old_stps[i])

        # multiply by initial Hessian
        # r/d is the final direction
        r = torch.mul(q, H_diag)
        for i in range(num_old):
            beta = old_stps[i].dot(r) * rho[i]
            r.add_(alpha[i] - beta, old_dirs[i])

        return r

    def curvature_update(self, flat_grad, eps=1e-2, damping=False):
        """
        Performs curvature update.

        Inputs:
            flat_grad (tensor): 1-D tensor of flattened gradient for computing
                gradient difference with previously stored gradient
            eps (float): constant for curvature pair rejection or damping (default: 1e-2)
            damping (bool): flag for using Powell damping (default: False)
        """

        assert len(self.param_groups) == 1

        # load parameters
        if eps <= 0:
            raise (ValueError("Invalid eps; must be positive."))

        group = self.param_groups[0]
        history_size = group["history_size"]
        debug = group["debug"]

        # variables cached in state (for tracing)
        state = self.state["global_state"]
        fail = state.get("fail")

        # check if line search failed
        if not fail:

            d = state.get("d")
            t = state.get("t")
            old_dirs = state.get("old_dirs")
            old_stps = state.get("old_stps")
            H_diag = state.get("H_diag")
            prev_flat_grad = state.get("prev_flat_grad")
            Bs = state.get("Bs")

            # compute y's
            y = flat_grad.sub(prev_flat_grad)
            s = d.mul(t)
            sBs = s.dot(Bs)
            ys = y.dot(s)  # y*s

            # update L-BFGS matrix
            if ys > eps * sBs or damping == True:

                # perform Powell damping
                if damping == True and ys < eps * sBs:
                    if debug:
                        print("Applying Powell damping...")
                    theta = ((1 - eps) * sBs) / (sBs - ys)
                    y = theta * y + (1 - theta) * Bs

                # updating memory
                if len(old_dirs) == history_size:
                    # shift history by one (limited-memory)
                    old_dirs.pop(0)
                    old_stps.pop(0)

                # store new direction/step
                old_dirs.append(s)
                old_stps.append(y)

                # update scale of initial Hessian approximation
                H_diag = ys / y.dot(y)  # (y*y)

                state["old_dirs"] = old_dirs
                state["old_stps"] = old_stps
                state["H_diag"] = H_diag

            else:
                # save skip
                state["curv_skips"] += 1
                if debug:
                    print("Curvature pair skipped due to failed criterion")

        else:
            # save skip
            state["fail_skips"] += 1
            if debug:
                print("Line search failed; curvature pair update skipped")

        return

    def _step(self, p_k, g_Ok, g_Sk=None, options=None):
        """
        Performs a single optimization step.

        Inputs:
            p_k (tensor): 1-D tensor specifying search direction
            g_Ok (tensor): 1-D tensor of flattened gradient over overlap O_k used
                            for gradient differencing in curvature pair update
            g_Sk (tensor): 1-D tensor of flattened gradient over full sample S_k
                            used for curvature pair damping or rejection criterion,
                            if None, will use g_Ok (default: None)
            options (dict): contains options for performing line search (default: None)

        Options for Armijo backtracking line search:
            'closure' (callable): reevaluates model and returns function value
            'current_loss' (tensor): objective value at current iterate (default: F(x_k))
            'gtd' (tensor): inner product g_Ok'd in line search (default: g_Ok'd)
            'eta' (tensor): factor for decreasing steplength > 0 (default: 2)
            'c1' (tensor): sufficient decrease constant in (0, 1) (default: 1e-4)
            'max_ls' (int): maximum number of line search steps permitted (default: 10)
            'interpolate' (bool): flag for using interpolation (default: True)
            'inplace' (bool): flag for inplace operations (default: True)
            'ls_debug' (bool): debugging mode for line search

        Options for Wolfe line search:
            'closure' (callable): reevaluates model and returns function value
            'current_loss' (tensor): objective value at current iterate (default: F(x_k))
            'gtd' (tensor): inner product g_Ok'd in line search (default: g_Ok'd)
            'eta' (float): factor for extrapolation (default: 2)
            'c1' (float): sufficient decrease constant in (0, 1) (default: 1e-4)
            'c2' (float): curvature condition constant in (0, 1) (default: 0.9)
            'max_ls' (int): maximum number of line search steps permitted (default: 10)
            'interpolate' (bool): flag for using interpolation (default: True)
            'inplace' (bool): flag for inplace operations (default: True)
            'ls_debug' (bool): debugging mode for line search

        Outputs (depends on line search):
          . No line search:
                t (float): steplength
          . Armijo backtracking line search:
                F_new (tensor): loss function at new iterate
                t (tensor): final steplength
                ls_step (int): number of backtracks
                closure_eval (int): number of closure evaluations
                desc_dir (bool): descent direction flag
                    True: p_k is descent direction with respect to the line search
                    function
                    False: p_k is not a descent direction with respect to the line
                    search function
                fail (bool): failure flag
                    True: line search reached maximum number of iterations, failed
                    False: line search succeeded
          . Wolfe line search:
                F_new (tensor): loss function at new iterate
                g_new (tensor): gradient at new iterate
                t (float): final steplength
                ls_step (int): number of backtracks
                closure_eval (int): number of closure evaluations
                grad_eval (int): number of gradient evaluations
                desc_dir (bool): descent direction flag
                    True: p_k is descent direction with respect to the line search
                    function
                    False: p_k is not a descent direction with respect to the line
                    search function
                fail (bool): failure flag
                    True: line search reached maximum number of iterations, failed
                    False: line search succeeded

        Notes:
          . If encountering line search failure in the deterministic setting, one
            should try increasing the maximum number of line search steps max_ls.

        """

        if options is None:
            options = {}
        assert len(self.param_groups) == 1

        # load parameter options
        group = self.param_groups[0]
        lr = group["lr"]
        line_search = group["line_search"]
        dtype = group["dtype"]
        debug = group["debug"]

        # variables cached in state (for tracing)
        state = self.state["global_state"]
        d = state.get("d")
        t = state.get("t")
        prev_flat_grad = state.get("prev_flat_grad")
        Bs = state.get("Bs")

        # keep track of nb of iterations
        state["n_iter"] += 1

        # set search direction
        d = p_k

        # modify previous gradient
        if prev_flat_grad is None:
            prev_flat_grad = g_Ok.clone()
        else:
            prev_flat_grad.copy_(g_Ok)

        # set initial step size
        t = lr

        # closure evaluation counter
        closure_eval = 0

        if g_Sk is None:
            g_Sk = g_Ok.clone()

        # perform Armijo backtracking line search
        if line_search == "Armijo":

            # load options
            if options:
                if "closure" not in options.keys():
                    raise (ValueError("closure option not specified."))
                else:
                    closure = options["closure"]

                if "gtd" not in options.keys():
                    gtd = g_Sk.dot(d)
                else:
                    gtd = options["gtd"]

                if "current_loss" not in options.keys():
                    F_k = closure()
                    closure_eval += 1
                else:
                    F_k = options["current_loss"]

                if "eta" not in options.keys():
                    eta = 2
                elif options["eta"] <= 0:
                    raise (ValueError("Invalid eta; must be positive."))
                else:
                    eta = options["eta"]

                if "c1" not in options.keys():
                    c1 = 1e-4
                elif options["c1"] >= 1 or options["c1"] <= 0:
                    raise (ValueError("Invalid c1; must be strictly between 0 and 1."))
                else:
                    c1 = options["c1"]

                if "max_ls" not in options.keys():
                    max_ls = 10
                elif options["max_ls"] <= 0:
                    raise (ValueError("Invalid max_ls; must be positive."))
                else:
                    max_ls = options["max_ls"]

                if "interpolate" not in options.keys():
                    interpolate = True
                else:
                    interpolate = options["interpolate"]

                if "inplace" not in options.keys():
                    inplace = True
                else:
                    inplace = options["inplace"]

                if "ls_debug" not in options.keys():
                    ls_debug = False
                else:
                    ls_debug = options["ls_debug"]

            else:
                raise (
                    ValueError(
                        "Options are not specified; need closure evaluating function."
                    )
                )

            # initialize values
            if interpolate:
                if torch.cuda.is_available():
                    F_prev = torch.tensor(np.nan, dtype=dtype).cuda()
                else:
                    F_prev = torch.tensor(np.nan, dtype=dtype)

            ls_step = 0
            t_prev = 0  # old steplength
            fail = False  # failure flag

            # begin print for debug mode
            if ls_debug:
                print(
                    "==================================== Begin Armijo line search ==================================="
                )
                print("F(x): %.8e  g*d: %.8e" % (F_k, gtd))

            # check if search direction is descent direction
            if gtd >= 0:
                desc_dir = False
                if debug:
                    print("Not a descent direction!")
            else:
                desc_dir = True

            # store values if not in-place
            if not inplace:
                current_params = self._copy_params()

            # update and evaluate at new point
            self._add_update(t, d)
            F_new = closure()
            closure_eval += 1

            # print info if debugging
            if ls_debug:
                print(
                    "LS Step: %d  t: %.8e  F(x+td): %.8e  F-c1*t*g*d: %.8e  F(x): %.8e"
                    % (ls_step, t, F_new, F_k + c1 * t * gtd, F_k)
                )

            # check Armijo condition
            while F_new > F_k + c1 * t * gtd or not is_legal(F_new):

                # check if maximum number of iterations reached
                if ls_step >= max_ls:
                    if inplace:
                        self._add_update(-t, d)
                    else:
                        self._load_params(current_params)

                    t = 0
                    F_new = closure()
                    closure_eval += 1
                    fail = True
                    break

                else:
                    # store current steplength
                    t_new = t

                    # compute new steplength

                    # if first step or not interpolating, then multiply by factor
                    if ls_step == 0 or not interpolate or not is_legal(F_new):
                        t = t / eta

                    # if second step, use function value at new point along with
                    # gradient and function at current iterate
                    elif ls_step == 1 or not is_legal(F_prev):
                        t = polyinterp(
                            np.array(
                                [
                                    [0, F_k.item(), gtd.item()],
                                    [t_new, F_new.item(), np.nan],
                                ]
                            )
                        )

                    # otherwise, use function values at new point, previous point,
                    # and gradient and function at current iterate
                    else:
                        t = polyinterp(
                            np.array(
                                [
                                    [0, F_k.item(), gtd.item()],
                                    [t_new, F_new.item(), np.nan],
                                    [t_prev, F_prev.item(), np.nan],
                                ]
                            )
                        )

                    # if values are too extreme, adjust t
                    if interpolate:
                        if t < 1e-3 * t_new:
                            t = 1e-3 * t_new
                        elif t > 0.6 * t_new:
                            t = 0.6 * t_new

                        # store old point
                        F_prev = F_new
                        t_prev = t_new

                    # update iterate and reevaluate
                    if inplace:
                        self._add_update(t - t_new, d)
                    else:
                        self._load_params(current_params)
                        self._add_update(t, d)

                    F_new = closure()
                    closure_eval += 1
                    ls_step += 1  # iterate

                    # print info if debugging
                    if ls_debug:
                        print(
                            "LS Step: %d  t: %.8e  F(x+td):   %.8e  F-c1*t*g*d: %.8e  F(x): %.8e"
                            % (ls_step, t, F_new, F_k + c1 * t * gtd, F_k)
                        )

            # store Bs
            if Bs is None:
                Bs = (g_Sk.mul(-t)).clone()
            else:
                Bs.copy_(g_Sk.mul(-t))

            # print final steplength
            if ls_debug:
                print("Final Steplength:", t)
                print(
                    "===================================== End Armijo line search ===================================="
                )

            state["d"] = d
            state["prev_flat_grad"] = prev_flat_grad
            state["t"] = t
            state["Bs"] = Bs
            state["fail"] = fail

            return F_new, t, ls_step, closure_eval, desc_dir, fail

        # perform weak Wolfe line search
        elif line_search == "Wolfe":

            # load options
            if options:
                if "closure" not in options.keys():
                    raise (ValueError("closure option not specified."))
                else:
                    closure = options["closure"]

                if "current_loss" not in options.keys():
                    F_k = closure()
                    closure_eval += 1
                else:
                    F_k = options["current_loss"]

                if "gtd" not in options.keys():
                    gtd = g_Sk.dot(d)
                else:
                    gtd = options["gtd"]

                if "eta" not in options.keys():
                    eta = 2
                elif options["eta"] <= 1:
                    raise (ValueError("Invalid eta; must be greater than 1."))
                else:
                    eta = options["eta"]

                if "c1" not in options.keys():
                    c1 = 1e-4
                elif options["c1"] >= 1 or options["c1"] <= 0:
                    raise (ValueError("Invalid c1; must be strictly between 0 and 1."))
                else:
                    c1 = options["c1"]

                if "c2" not in options.keys():
                    c2 = 0.9
                elif options["c2"] >= 1 or options["c2"] <= 0:
                    raise (ValueError("Invalid c2; must be strictly between 0 and 1."))
                elif options["c2"] <= c1:
                    raise (ValueError("Invalid c2; must be strictly larger than c1."))
                else:
                    c2 = options["c2"]

                if "max_ls" not in options.keys():
                    max_ls = 10
                elif options["max_ls"] <= 0:
                    raise (ValueError("Invalid max_ls; must be positive."))
                else:
                    max_ls = options["max_ls"]

                if "interpolate" not in options.keys():
                    interpolate = True
                else:
                    interpolate = options["interpolate"]

                if "inplace" not in options.keys():
                    inplace = True
                else:
                    inplace = options["inplace"]

                if "ls_debug" not in options.keys():
                    ls_debug = False
                else:
                    ls_debug = options["ls_debug"]

            else:
                raise (
                    ValueError(
                        "Options are not specified; need closure evaluating function."
                    )
                )

            # initialize counters
            ls_step = 0
            grad_eval = 0  # tracks gradient evaluations
            t_prev = 0  # old steplength

            # initialize bracketing variables and flag
            alpha = 0
            beta = float("Inf")
            fail = False

            # initialize values for line search
            if interpolate:
                F_a = F_k
                g_a = gtd

                if torch.cuda.is_available():
                    F_b = torch.tensor(np.nan, dtype=dtype).cuda()
                    g_b = torch.tensor(np.nan, dtype=dtype).cuda()
                else:
                    F_b = torch.tensor(np.nan, dtype=dtype)
                    g_b = torch.tensor(np.nan, dtype=dtype)

            # begin print for debug mode
            if ls_debug:
                print(
                    "==================================== Begin Wolfe line search ===================================="
                )
                print("F(x): %.8e  g*d: %.8e" % (F_k, gtd))

            # check if search direction is descent direction
            if gtd >= 0:
                desc_dir = False
                if debug:
                    print("Not a descent direction!")
            else:
                desc_dir = True

            # store values if not in-place
            if not inplace:
                current_params = self._copy_params()

            # update and evaluate at new point
            self._add_update(t, d)
            F_new = closure()
            closure_eval += 1

            # main loop
            while True:

                # check if maximum number of line search steps have been reached
                if ls_step >= max_ls:
                    if inplace:
                        self._add_update(-t, d)
                    else:
                        self._load_params(current_params)

                    t = 0
                    F_new = closure()
                    F_new.backward()
                    g_new = self._gather_flat_grad()
                    closure_eval += 1
                    grad_eval += 1
                    fail = True
                    break

                # print info if debugging
                if ls_debug:
                    print(
                        "LS Step: %d  t: %.8e  alpha: %.8e  beta: %.8e"
                        % (ls_step, t, alpha, beta)
                    )
                    print(
                        "Armijo:  F(x+td): %.8e  F-c1*t*g*d: %.8e  F(x): %.8e"
                        % (F_new, F_k + c1 * t * gtd, F_k)
                    )

                # check Armijo condition
                if F_new > F_k + c1 * t * gtd:

                    # set upper bound
                    beta = t
                    t_prev = t

                    # update interpolation quantities
                    if interpolate:
                        F_b = F_new
                        if torch.cuda.is_available():
                            g_b = torch.tensor(np.nan, dtype=dtype).cuda()
                        else:
                            g_b = torch.tensor(np.nan, dtype=dtype)

                else:

                    # compute gradient
                    F_new.backward()
                    g_new = self._gather_flat_grad()
                    grad_eval += 1
                    gtd_new = g_new.dot(d)

                    # print info if debugging
                    if ls_debug:
                        print(
                            "Wolfe: g(x+td)*d: %.8e  c2*g*d: %.8e  gtd: %.8e"
                            % (gtd_new, c2 * gtd, gtd)
                        )

                    # check curvature condition
                    if gtd_new < c2 * gtd:

                        # set lower bound
                        alpha = t
                        t_prev = t

                        # update interpolation quantities
                        if interpolate:
                            F_a = F_new
                            g_a = gtd_new

                    else:
                        break

                # compute new steplength

                # if first step or not interpolating, then bisect or multiply by factor
                if not interpolate or not is_legal(F_b):
                    if beta == float("Inf"):
                        t = eta * t
                    else:
                        t = (alpha + beta) / 2.0

                # otherwise interpolate between a and b
                else:
                    t = polyinterp(
                        np.array(
                            [
                                [alpha, F_a.item(), g_a.item()],
                                [beta, F_b.item(), g_b.item()],
                            ]
                        )
                    )

                    # if values are too extreme, adjust t
                    if beta == float("Inf"):
                        if t > 2 * eta * t_prev:
                            t = 2 * eta * t_prev
                        elif t < eta * t_prev:
                            t = eta * t_prev
                    else:
                        if t < alpha + 0.2 * (beta - alpha):
                            t = alpha + 0.2 * (beta - alpha)
                        elif t > (beta - alpha) / 2.0:
                            t = (beta - alpha) / 2.0

                    # if we obtain nonsensical value from interpolation
                    if t <= 0:
                        t = (beta - alpha) / 2.0

                # update parameters
                if inplace:
                    self._add_update(t - t_prev, d)
                else:
                    self._load_params(current_params)
                    self._add_update(t, d)

                # evaluate closure
                F_new = closure()
                closure_eval += 1
                ls_step += 1

            # store Bs
            if Bs is None:
                Bs = (g_Sk.mul(-t)).clone()
            else:
                Bs.copy_(g_Sk.mul(-t))

            # print final steplength
            if ls_debug:
                print("Final Steplength:", t)
                print(
                    "===================================== End Wolfe line search ====================================="
                )

            state["d"] = d
            state["prev_flat_grad"] = prev_flat_grad
            state["t"] = t
            state["Bs"] = Bs
            state["fail"] = fail

            return F_new, g_new, t, ls_step, closure_eval, grad_eval, desc_dir, fail

        else:

            # perform update
            self._add_update(t, d)

            # store Bs
            if Bs is None:
                Bs = (g_Sk.mul(-t)).clone()
            else:
                Bs.copy_(g_Sk.mul(-t))

            state["d"] = d
            state["prev_flat_grad"] = prev_flat_grad
            state["t"] = t
            state["Bs"] = Bs
            state["fail"] = False

            return t

    def step(self, p_k, g_Ok, g_Sk=None, options={}):
        return self._step(p_k, g_Ok, g_Sk, options)


class FullBatchLBFGS(LBFGS):
    """
    Implements full-batch or deterministic L-BFGS algorithm. Compatible with
    Powell damping. Can be used when evaluating a deterministic function and
    gradient. Wraps the LBFGS optimizer. Performs the two-loop recursion,
    updating, and curvature updating in a single step.

    Implemented by: Hao-Jun Michael Shi and Dheevatsa Mudigere
    Last edited 11/15/18.

    Warnings:
      . Does not support per-parameter options and parameter groups.
      . All parameters have to be on a single device.

    Inputs:
        lr (float): steplength or learning rate (default: 1)
        history_size (int): update history size (default: 10)
        line_search (str): designates line search to use (default: 'Wolfe')
            Options:
                'None': uses steplength designated in algorithm
                'Armijo': uses Armijo backtracking line search
                'Wolfe': uses Armijo-Wolfe bracketing line search
        dtype: data type (default: torch.float)
        debug (bool): debugging mode

    """

    def __init__(
        self,
        params,
        lr=1,
        history_size=10,
        line_search="Wolfe",
        dtype=torch.float,
        debug=False,
    ):
        super(FullBatchLBFGS, self).__init__(
            params, lr, history_size, line_search, dtype, debug
        )

    def step(self, options=None):
        """
        Performs a single optimization step.

        Inputs:
            options (dict): contains options for performing line search (default: None)

        General Options:
            'eps' (float): constant for curvature pair rejection or damping (default: 1e-2)
            'damping' (bool): flag for using Powell damping (default: False)

        Options for Armijo backtracking line search:
            'closure' (callable): reevaluates model and returns function value
            'current_loss' (tensor): objective value at current iterate (default: F(x_k))
            'gtd' (tensor): inner product g_Ok'd in line search (default: g_Ok'd)
            'eta' (tensor): factor for decreasing steplength > 0 (default: 2)
            'c1' (tensor): sufficient decrease constant in (0, 1) (default: 1e-4)
            'max_ls' (int): maximum number of line search steps permitted (default: 10)
            'interpolate' (bool): flag for using interpolation (default: True)
            'inplace' (bool): flag for inplace operations (default: True)
            'ls_debug' (bool): debugging mode for line search

        Options for Wolfe line search:
            'closure' (callable): reevaluates model and returns function value
            'current_loss' (tensor): objective value at current iterate (default: F(x_k))
            'gtd' (tensor): inner product g_Ok'd in line search (default: g_Ok'd)
            'eta' (float): factor for extrapolation (default: 2)
            'c1' (float): sufficient decrease constant in (0, 1) (default: 1e-4)
            'c2' (float): curvature condition constant in (0, 1) (default: 0.9)
            'max_ls' (int): maximum number of line search steps permitted (default: 10)
            'interpolate' (bool): flag for using interpolation (default: True)
            'inplace' (bool): flag for inplace operations (default: True)
            'ls_debug' (bool): debugging mode for line search

        Outputs (depends on line search):
          . No line search:
                t (float): steplength
          . Armijo backtracking line search:
                F_new (tensor): loss function at new iterate
                t (tensor): final steplength
                ls_step (int): number of backtracks
                closure_eval (int): number of closure evaluations
                desc_dir (bool): descent direction flag
                    True: p_k is descent direction with respect to the line search
                    function
                    False: p_k is not a descent direction with respect to the line
                    search function
                fail (bool): failure flag
                    True: line search reached maximum number of iterations, failed
                    False: line search succeeded
          . Wolfe line search:
                F_new (tensor): loss function at new iterate
                g_new (tensor): gradient at new iterate
                t (float): final steplength
                ls_step (int): number of backtracks
                closure_eval (int): number of closure evaluations
                grad_eval (int): number of gradient evaluations
                desc_dir (bool): descent direction flag
                    True: p_k is descent direction with respect to the line search
                    function
                    False: p_k is not a descent direction with respect to the line
                    search function
                fail (bool): failure flag
                    True: line search reached maximum number of iterations, failed
                    False: line search succeeded

        Notes:
          . If encountering line search failure in the deterministic setting, one
            should try increasing the maximum number of line search steps max_ls.

        """

        # load options for damping and eps
        if "damping" not in options.keys():
            damping = False
        else:
            damping = options["damping"]

        if "eps" not in options.keys():
            eps = 1e-2
        else:
            eps = options["eps"]

        # gather gradient
        grad = self._gather_flat_grad()
        # /print(grad)
        # exit()
        # update curvature if after 1st iteration
        state = self.state["global_state"]
        if state["n_iter"] > 0:
            self.curvature_update(grad, eps, damping)

        # compute search direction
        p = self.two_loop_recursion(-grad)

        # take step
        return self._step(p, grad, options=options)

# execute script
if __name__ == '__main__':
    main()