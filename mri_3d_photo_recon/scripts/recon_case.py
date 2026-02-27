# Imports
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from photo_reconstruction import utils

# Constants we shouldn't have to touch
RESOLUTIONS = [4, 2, 1]  # Resolutions at which we work (in mm)
# Solves problems with machine learning imputation
# RESOLUTION_OUTPUT = 0.5  # Resolution for final 3D recon (in mm); will be ignored (set to 1mm) if weights are provided
RESOLUTION_OUTPUT = 1.0
STEPS = [250, 150, 75]  # Optimization steps at every resolution (per pass)

TOL = 1e-6  # Tolerance / threshold for convergence
PAD_AP = 3  # Padding in A-P axis ensures that registered 3D mask remains 100% in FOV
THRESHOLD_FG = 5.0 # for ML normalization
UNSHARP_SIGMA = 1.0
UNSHARP_AMOUNT = 1.0
LR = 0.02
# OPTIMIZER_TYPE = 'Adam'
OPTIMIZER_TYPE = 'NAdam'

# ================================================================================================
#                                         Main Entrypoint
# ================================================================================================
def main():

    # Get arguments from command line and decide what transforms to allow and with what regularizers, control
    # point spacings, etc. We also play a little trick to help the first/last few slices remain in the mesh
    arguments, allow_z_stretch, allow_affine_mri = utils.adjust_settings(utils.get_arguments())
    # arguments.slice_thickness = (0.8*arguments.slice_thickness) if allow_z_stretch else arguments.slice_thickness

    # only import packages if arguments are correct
    from nibabel.freesurfer import write_geometry
    import numpy as np
    from scipy.ndimage import gaussian_filter, binary_fill_holes
    import torch
    from photo_reconstruction import image_utils, mesh_utils
    from photo_reconstruction.photo_aligner import photo_aligner
    from photo_reconstruction.optimization import FullBatchLBFGS
    from datetime import datetime
    import random

    # Keep track of time
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Current Time =", current_time)

    # Sort out gpu and threads
    device = utils.configure_gpu_and_cpu(arguments.gpu, arguments.threads)

    # Reproducibility: https://discuss.pytorch.org/t/reproducibility-with-all-the-bells-and-whistles/81097
    seed = 0
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.cuda.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    # Data types and seeds, right off the bat
    torch.set_default_dtype(torch.float32)
    utils.seed_all(seed)

    # Optimizer
    if OPTIMIZER_TYPE=='NAdam':
        optimizer_type = torch.optim.NAdam
    elif OPTIMIZER_TYPE=='Adam':
        optimizer_type = torch.optim.NAdam
    else:
        raise Exception('Optimizer not supported: ' + OPTIMIZER_TYPE)

    # Get prefix for FreeSurfer commands (sources freesurfer)
    fsprefix = utils.get_fs_prefix()

    # Make output directory if needed
    output_directory = arguments.output_directory + "/"
    if os.path.isdir(output_directory) is False:
        os.makedirs(output_directory, exist_ok=True)

    # Extract slices from photographs, crop to equal size, and make distance maps for masks
    d_i, d_s = utils.get_photo_and_seg_lists(arguments.input_photo_dir, arguments.input_segmentation_dir)
    Iorig, Morig = image_utils.read_images_and_masks(d_i, d_s, reverse_ap=arguments.posterior_to_anterior,
                                                     crop_margin=np.round(5 / arguments.photo_resolution),
                                                     ndilations=int(1.0 / arguments.photo_resolution),
                                                     equalize_images=arguments.equalize_images)
    areas = areas = np.array(Morig).sum(axis=(1, 2))
    for i in range(len(Morig)):
        Morig[i] = binary_fill_holes(Morig[i]>0).astype(Morig[i].dtype)
    Morig = image_utils.make_distance_transforms(Morig, arguments.photo_resolution)

    # get thicknesses from weights if needed, and estimate a-p coordinate shifts in mm
    # Crucial: first shift applies to first slab and last shift to slab N+1 (ie the first padding/empty slice after the last tissue slab)
    # (the alternative would be that the first shift applies to the empty slab right before the first tissue slab)
    if arguments.weights is None:
        print('No weight file provided; assuming constant thickness')
        y_shifts = None
    else:
        thicknesses = image_utils.estimate_thicknesses_from_weights(arguments.weights, areas, arguments.slice_thickness,
                                                                    arguments.thickness_cap, reverse_ap=arguments.posterior_to_anterior)
        y = arguments.slice_thickness * np.arange(PAD_AP)
        y = np.concatenate([y, y[-1] + thicknesses.cumsum()])
        y = np.concatenate([y, y[-1] + arguments.slice_thickness * np.arange(1, PAD_AP)])
        yref = arguments.slice_thickness * np.arange(len(y))
        y_shifts = y - yref

    # Resample to output resolution
    Iorig, Morig = image_utils.resample_images_and_masks(Iorig, Morig, arguments.photo_resolution, RESOLUTION_OUTPUT,
                                                         stretch_lr=arguments.initial_stretch_factor_lr_photos)

    # Build original 3D volumes with a bit of anterior-posterior padding.
    # Crucially: padded planes in distance maps must have high values
    Iorig = np.pad(np.stack(Iorig, axis=2), ((0,0),(0,0),(PAD_AP,PAD_AP),(0,0)))
    Morig = np.pad(np.stack(Morig, axis=2), ((0,0),(0,0),(PAD_AP,PAD_AP)))
    for k in range(PAD_AP):
        Morig[:, :, PAD_AP-k-1] =  - (k+1) * arguments.slice_thickness * np.ones_like(Morig[:, :, PAD_AP-k])
        Morig[:, :, -PAD_AP+k] = - (k+1) * arguments.slice_thickness * np.ones_like(Morig[:, :, -PAD_AP+k-1])
    siz = Iorig.shape[0:2]

    # Build resolution pyramid and center in origin
    Is, Ms, Affs, aff_orig = image_utils.make_pyramid(Iorig, Morig, RESOLUTIONS, RESOLUTION_OUTPUT, arguments.slice_thickness, arguments.posterior_side)
    aff_orig, cog_photo_ras = image_utils.center_affine_matrix_in_origin(Morig>0.5, aff_orig)
    for s in range(len(Affs)):
        Affs[s][:-1, -1] = Affs[s][:-1, -1] - np.squeeze(cog_photo_ras)

    # Apply masks now that all resampling is done
    Iorig[Morig<0] = 0
    for s in range(len(Is)):
        Is[s][Ms[s]<0] = 0

    # Preprocess MRI and segmentation, if available
    if (arguments.ref_mri is None) and (arguments.ref_mri_synthsr is None):
        if (arguments.ref_mesh is None):
            print('Neither 3D MRI nor mesh provided to reconstruct; we will use MNI as reference')
            REF, REF_SSEG, REFaff = utils.prepare_reference_volumes(os.path.join(os.path.dirname(__file__), '../mni/mni.nii.gz'), None, None,
                                                                    os.path.join(os.path.dirname(__file__), '../mni/mni.seg.nii.gz'),
                                                                    output_directory, fsprefix, arguments.threads)
        else:
            print('3D MRI reference not provided')
            REF = REF_SSEG = REFmask = REFaff = cog_mri_ras = None
    else:
        REF, REF_SSEG, REFaff = utils.prepare_reference_volumes(arguments.ref_mri, arguments.ref_mri_synthsr, arguments.low_field_synthsr,
                                                                arguments.ref_mri_synthseg, output_directory, fsprefix, arguments.threads)
    if REF is not None:
        REFmask = image_utils.compute_mask_from_synthseg(REF_SSEG, arguments.hemisphere)
        REFaff, cog_mri_ras = image_utils.center_affine_matrix_in_origin(REFmask, REFaff)
        
    # Read and preprocess mesh, if available.
    # We play a little trick to make sure that the first/last few slices do not get ignored,
    # by adding points randomly around the cloud with lower weight (improves capture range)
    # Note to self: I disabled it as it could stretch the distance transform / edge in a funny way
    if (arguments.ref_mesh is None):
        Pmesh = Dmesh = TRImesh = meta_mesh = None
    else:
        Pmesh, TRImesh, meta_mesh = mesh_utils.read_and_reorient_mesh(arguments.ref_mesh, arguments.mesh_reorient_with_indices, fsprefix, output_directory, arguments.hemisphere)
        Pmesh[:,0] = arguments.stretch_factor_lr_mesh * Pmesh[:,0]
        nv_orig = Pmesh.shape[0]
        Dmesh = np.zeros(nv_orig)
        print('Computing approximate distance function for mesh')
        margin = 15
        mini = (Pmesh.min(axis=0) - margin).astype(np.int32)
        maxi = (Pmesh.max(axis=0) + margin + 1).astype(np.int32)
        np_target = 50000
        # I used to make delta a function of the number of mesh vertices, but this led to too few points sometimes
        delta = 1.0 
        skip = int(np.floor(Pmesh.shape[0]/np_target))
        points = torch.tensor(Pmesh[::skip,:], device=device, dtype=torch.float)
        xs = torch.arange(mini[0], maxi[0], delta, device=device, dtype=torch.float32)
        ys = torch.arange(mini[1], maxi[1], delta, device=device, dtype=torch.float32)
        zs = torch.arange(mini[2], maxi[2], delta, device=device, dtype=torch.float32)
        grid = torch.stack(torch.meshgrid(xs, ys, zs, indexing="ij"), -1).reshape(-1, 3)
        df_vals = torch.empty(grid.shape[0], device=device)
        batch_size = 10000
        for i in range(0, grid.shape[0], batch_size):
            q = grid[i:i + batch_size]  # (B,3)
            dists = torch.cdist(q, points)  # (B,N)
            df_vals[i:i + batch_size] = dists.min(dim=1).values
        Pmesh = np.concatenate([Pmesh, grid.detach().cpu().numpy()], axis=0)
        Dmesh = np.concatenate([Dmesh, df_vals.detach().cpu().numpy()], axis=0)
        print('Done computing distances')
        del points, grid, df_vals

    ########################################################

    print("Optimization")

    # Initialize
    t = theta = shear = scaling = sz = t_mri = theta_mri = shear_mri = s_mri = t_mesh = theta_mesh = field = field3d = None

    # Go over resolutions / modes
    print("We will be running 3 modes: rigid, affine, and nonlinear")

    for mode_idx in range(3):
        if mode_idx == 0:
            print("\n###    First pass: no scaling / shearing allowed     ###")
            allow_scaling_and_shear = allow_nonlin = allow_nonlin_mri = alternate_optimization = False
        elif mode_idx==1:
            print("\n###    Second pass: scaling / shearing is allowed    ###")
            allow_scaling_and_shear = True
            allow_nonlin = allow_nonlin_mri = alternate_optimization = False
        else:
            print("\n###    Third pass: nonlinear deformation is allowed ###")
            allow_scaling_and_shear = allow_nonlin = allow_nonlin_mri = True
            alternate_optimization = False

        learning_rate_adam = LR
        trigger_limit = 15 if alternate_optimization else 5
        for res in range(len(RESOLUTIONS)):

            print("\nWorking on resolution %d of %d (%.2f mm): %d iterations "
                % (res + 1, len(RESOLUTIONS), RESOLUTIONS[res], STEPS[res]) )

            # Smooth MRI and mask to target resolution
            if REF is None:
                REFsmooth = REFmaskSmooth = None
            else:
                volres = np.sqrt(np.sum(REFaff[:, :-1] ** 2, axis=0))
                sigmas = 0.5 * RESOLUTIONS[res] / volres
                REFsmooth = gaussian_filter(REF, sigmas)
                REFmaskSmooth = gaussian_filter(REFmask.astype(np.float32), sigmas)

            # Build aligner model
            model = photo_aligner(
                Is[res],
                Ms[res],
                Affs[res],
                REFsmooth,
                REFmaskSmooth,
                REFaff,
                Pmesh,
                Dmesh,
                pixel_size=RESOLUTIONS[res],
                t_ini=t,
                theta_ini=theta,
                shear_ini=shear,
                scaling_ini=scaling,
                sz_ini=sz,
                y_shifts=y_shifts,
                field_ini=field,
                t_mri_ini=t_mri,
                theta_mri_ini=theta_mri,
                shear_mri_ini=shear_mri,
                s_mri_ini=s_mri,
                field3d_ini=field3d,
                t_mesh_ini=t_mesh,
                theta_mesh_ini=theta_mesh,
                allow_scaling_and_shear=allow_scaling_and_shear,
                allow_sz=allow_z_stretch,
                allow_nonlin=allow_nonlin,
                allow_affine_mri = allow_affine_mri,
                allow_nonlin_mri=allow_nonlin_mri,
                cp_spacing2d=arguments.cp_spacing_2d,
                cp_spacing3d=arguments.cp_spacing_3d,
                k_lncc_mri=arguments.k_lncc_mri,
                k_dice_mri=arguments.k_dice_mri,
                k_dif_slice_loss=arguments.k_dif_slice_loss,
                k_mesh_loss=arguments.k_mesh_loss,
                k_regularizer=arguments.k_regularizer,
                k_regularizer_sz=arguments.k_regularizer_sz,
                k_regularizer_nonlin=arguments.k_regularizer_nonlin,
                k_regularizer_nonlin3d=arguments.k_regularizer_nonlin3d,
                pad_ignore=PAD_AP,
                device=device
            )

            loss_old = 1e10
            for optimizer_mode in range(2):
                # Optimization with (N)Adam
                if optimizer_mode==0:
                    nsteps = STEPS[res]
                    if alternate_optimization: # doesn't make much sense with Adam, I think, but I kept it there
                        optimizer_photos = optimizer_type(model.parameters_photos(), lr=learning_rate_adam)
                        optimizer_mri = None if (REF is None) else optimizer_type(model.parameters_mri(), lr=learning_rate_adam)
                        optimizer_mesh = None if (Pmesh is None) else optimizer_type(model.parameters_mesh(), lr=learning_rate_adam)
                    else:
                        optimizer = optimizer_type(model.parameters(), lr=learning_rate_adam)
                    options = options_photos = options_mri = options_mesh = None
                elif arguments.skip_bfgs:
                    # Do nothing
                    nsteps = 0
                else:
                    # Finetuning with LBFGS
                    nsteps = 15 # good for alternate mode
                    print('\n   We fineture with 15 LBFGS iterations')
                    if alternate_optimization:
                        optimizer_photos = FullBatchLBFGS(model.parameters_photos())
                        optimizer_mri = None if (REF is None) else FullBatchLBFGS(model.parameters_mri())
                        optimizer_mesh = None if (Pmesh is None) else FullBatchLBFGS(model.parameters_mesh())
                        def closure_photos():
                            optimizer_photos.zero_grad(); loss = model()[0]; return loss
                        def closure_mri():
                            optimizer_mri.zero_grad(); loss = model()[0]; return loss
                        def closure_mesh():
                            optimizer_mesh.zero_grad(); loss = model()[0]; return loss
                        options_photos = {"closure": closure_photos, "current_loss": loss, "max_ls": 75}
                        options_mri = {"closure": closure_mri, "current_loss": loss, "max_ls": 75}
                        options_mesh = {"closure": closure_mesh, "current_loss": loss, "max_ls": 75}
                    else:
                        optimizer = FullBatchLBFGS(model.parameters())
                        def closure():
                            optimizer.zero_grad(); loss = model()[0]; return loss
                        options = {"closure": closure, "current_loss": loss, "max_ls": 75}

                # Iterate
                trigger_times = 0
                for epoch in range(nsteps):
                    # Compute loss with forward pass
                    loss = model()[0]
                    # Alternating scheme depends on what modalities are available
                    if alternate_optimization:
                        if (REF is not None) and (Pmesh is not None):
                            if epoch % 15 < 5:
                                optimizer_photos.zero_grad(); loss.backward(); optimizer_photos.step(options_photos)
                            elif epoch % 15 < 10:
                                optimizer_mri.zero_grad(); loss.backward(); optimizer_mri.step(options_mri)
                            else:
                                optimizer_mesh.zero_grad(); loss.backward(); optimizer_mesh.step(options_mesh)
                        else:
                            if epoch % 10 < 5:
                                optimizer_photos.zero_grad(); loss.backward(); optimizer_photos.step(options_photos)
                            else:
                                if (REF is None):
                                    optimizer_mesh.zero_grad(); loss.backward(); optimizer_mesh.step(options_mesh)
                                else:
                                    optimizer_mri.zero_grad(); loss.backward(); optimizer_mri.step(options_mri)
                    else:
                        optimizer.zero_grad(); loss.backward(); optimizer.step(options)

                    # print step info
                    loss_cpu = loss.cpu().detach().numpy()
                    print("   Step %d, loss = %.10f" % (epoch + 1, loss_cpu), flush=True, end='\r')

                    # Exit condition
                    if ((loss_old - loss_cpu) < TOL):
                        trigger_times += 1
                        if (trigger_times >= trigger_limit) and (optimizer_mode==0):
                            print("\n   Decrease in loss below tolerance limit for the last " +  str(trigger_limit) + " steps; dividing learning rate by 2")
                            learning_rate_adam *= 0.5
                            break
                    else:
                        trigger_times = 0
                    loss_old = loss_cpu

            # Retrieve model parameters
            print(' ')
            t = model.t.cpu().detach().numpy()
            theta = model.theta.cpu().detach().numpy()
            shear = model.shear.cpu().detach().numpy()
            scaling = model.scaling.cpu().detach().numpy()
            sz = model.sz.cpu().detach().numpy()
            field = None if (model.field is None) else model.field.cpu().detach().numpy()
            t_mri = None if (model.t_mri is None) else model.t_mri.cpu().detach().numpy()
            theta_mri = None if (model.theta_mri is None) else model.theta_mri.cpu().detach().numpy()
            shear_mri = None if (model.shear_mri is None) else model.shear_mri.cpu().detach().numpy()
            s_mri = None if (model.s_mri is None) else model.s_mri.cpu().detach().numpy()
            field3d = None if (model.field3d is None) else model.field3d.cpu().detach().numpy()
            t_mesh = None if (model.t_mesh is None) else model.t_mesh.cpu().detach().numpy()
            theta_mesh = None if (model.theta_mesh is None) else model.theta_mesh.cpu().detach().numpy()

            # In the last resolution level, retrieve results before deleting the model
            if res == (len(RESOLUTIONS) - 1) and mode_idx == 2:
                # We apply the deformations to the higher resolution photos
                model_hr = photo_aligner(
                    Iorig[...,0],
                    Morig,
                    aff_orig,
                    REF,
                    None if REF is None else REFmask.astype(np.float32),
                    REFaff,
                    Pmesh,
                    Dmesh,
                    pixel_size=RESOLUTION_OUTPUT,
                    t_ini=t,
                    theta_ini=theta,
                    shear_ini=shear,
                    scaling_ini=scaling,
                    sz_ini=sz,
                    y_shifts=y_shifts,
                    field_ini=field,
                    t_mri_ini=t_mri,
                    theta_mri_ini=theta_mri,
                    shear_mri_ini=shear_mri,
                    s_mri_ini=s_mri,
                    field3d_ini=field3d,
                    t_mesh_ini=t_mesh,
                    theta_mesh_ini=theta_mesh,
                    allow_scaling_and_shear=allow_scaling_and_shear,
                    allow_sz=allow_z_stretch,
                    allow_nonlin=allow_nonlin,
                    allow_affine_mri=allow_affine_mri,
                    allow_nonlin_mri=allow_nonlin_mri,
                    cp_spacing2d=arguments.cp_spacing_2d,
                    cp_spacing3d=arguments.cp_spacing_3d,
                    k_lncc_mri=arguments.k_lncc_mri,
                    k_dice_mri=arguments.k_dice_mri,
                    k_dif_slice_loss=arguments.k_dif_slice_loss,
                    k_mesh_loss=arguments.k_mesh_loss,
                    k_regularizer=arguments.k_regularizer,
                    k_regularizer_sz=arguments.k_regularizer_sz,
                    k_regularizer_nonlin=arguments.k_regularizer_nonlin,
                    k_regularizer_nonlin3d=arguments.k_regularizer_nonlin3d,
                    pad_ignore=PAD_AP,
                    device=device
                )
                (_, photo_resampled_r, photo_aff, mri_aff_combined, Rt, TvoxPhotos, mri_resampled, Tmesh, grids_new_mri_nonlin) = model_hr()
                model_hr.photo_vol = torch.Tensor(Iorig[:, :, :, 1].copy()).to(device)
                model_hr.photo_rearranged = torch.unsqueeze(torch.unsqueeze(model_hr.photo_vol, dim=0), dim=0).to(model_hr.device)
                aux = model_hr(); photo_resampled_g = aux[1]
                model_hr.photo_vol = torch.Tensor(Iorig[:, :, :, 2].copy()).to(device)
                model_hr.photo_rearranged = torch.unsqueeze(torch.unsqueeze(model_hr.photo_vol, dim=0), dim=0).to(model_hr.device)
                aux = model_hr(); photo_resampled_b = aux[1]
                photo_resampled = torch.stack([photo_resampled_r, photo_resampled_g, photo_resampled_b], dim=-1)

                # photo_resampled = photo_resampled.cpu().detach().numpy() # do not detach yet, we'll use later
                photo_aff = photo_aff.cpu().detach().numpy()
                Rt = None if Rt is None else Rt.cpu().detach().numpy()
                TvoxPhotos = None if TvoxPhotos is None else TvoxPhotos.cpu().detach().numpy()
                if REF is None:
                    mri_resampled = mri_aff_combined = None
                else:
                    mri_resampled = mri_resampled.cpu().detach().numpy()
                    mri_aff_combined = mri_aff_combined.cpu().detach().numpy()
                if Pmesh is None:
                    Tmesh = None
                else:
                    Tmesh = Tmesh.cpu().detach().numpy()
                if grids_new_mri_nonlin is not None:
                    grids_new_mri_nonlin = np.squeeze(grids_new_mri_nonlin.cpu().detach().numpy())
                    if y_shifts is None:
                        grids_new_mri_nonlin_no_shift = None
                    else:
                        model_hr.y_shifts = None
                        aux = model_hr()
                        grids_new_mri_nonlin_no_shift = np.squeeze(aux[8].cpu().detach().numpy())
                del model_hr

            # Free up memory
            if alternate_optimization:
                del optimizer_mri
                del optimizer_mesh
                del optimizer_photos
            else:
                del optimizer
            del model

    ########################################################

    # Write final outputs (run ML interpolation if uneven spacing!)
    cmd = 'freeview '
    if y_shifts is None:
        print('\nWriting uncorrected photo recon to disk')
        cmd += (' -v ' + output_directory + '/photo_recon.orig.mgz:rgb=1 ')
        image_utils.MRIwrite(photo_resampled.cpu().detach().numpy().clip(0,255), photo_aff, output_directory + '/photo_recon.orig.mgz', dtype=np.uint8)

    print('\nMachine learning interpolation')
    from photo_reconstruction.machine_learning_utils import photo_imputation
    LINEAR, PRED, affnew = photo_imputation(photo_resampled, photo_aff, y_shifts, sz, THRESHOLD_FG, UNSHARP_SIGMA, UNSHARP_AMOUNT,
                                    os.path.join(os.environ.get('FREESURFER_HOME_FSPYTHON'), 'models/photo_imputation_unet.pth'), device)
    image_utils.MRIwrite(LINEAR.cpu().detach().numpy().clip(0,255), affnew, output_directory + '/photo_recon.trilinear.mgz', dtype=np.uint8)
    image_utils.MRIwrite(PRED.cpu().detach().numpy().clip(0,255), affnew, output_directory + '/photo_recon.machine_learning.mgz', dtype=np.uint8)
    cmd += (' -v ' + output_directory + '/photo_recon.trilinear.mgz:rgb=1 ')
    cmd += (' -v ' + output_directory + '/photo_recon.machine_learning.mgz:rgb=1 ')

    print('\nWriting deformed references to disk')
    if REF is not None:
        if y_shifts is None:
            cmd += (' -v ' + output_directory + '/mri.deformed.photo_space.mgz ')
            image_utils.MRIwrite(mri_resampled, photo_aff, output_directory + '/mri.deformed.photo_space.mgz')
        REFaff[:-1, -1] += np.squeeze(cog_mri_ras)
        if y_shifts is None:
            RAS, RASres = image_utils_refine_with_seg.computeRAS(grids_new_mri_nonlin, REF.shape, REFaff, photo_aff, LINEAR.shape, affnew, device)
            image_utils.MRIwrite(RAS, photo_aff, output_directory + '/field.photo_space.mgz')
            image_utils.MRIwrite(RASres, affnew, output_directory + '/field.1mm.mgz')
        else:
            RAS, RASres = image_utils_refine_with_seg.computeRAS(grids_new_mri_nonlin_no_shift, REF.shape, REFaff, photo_aff, LINEAR.shape, affnew, device)
            image_utils.MRIwrite(RASres, affnew, output_directory + '/field.1mm.mgz')
        REFres, _ = image_utils.deform(REF, REFaff, RASres, device, mode='linear')
        REFmaskRes, _ = image_utils.deform(REFmask.astype(np.float32), REFaff, RASres, device, mode='linear')
        REFres *= REFmaskRes
        image_utils.MRIwrite(REFres, affnew, output_directory + '/mri.deformed.mgz')
        cmd += (' -v ' + output_directory + '/mri.deformed.mgz ')

        if arguments.input_roi_dir is not None:
            import glob
            g = sorted(glob.glob(arguments.input_roi_dir + '/*.nii')) + sorted(glob.glob(arguments.input_roi_dir + '/*.nii.gz')) + sorted(glob.glob(arguments.input_roi_dir + '/*.mgz'))
            print('Deforming ROIs')
            for i in range(len(g)):
                print('  ROI ' + str(i+1) + ' of ' + str(len(g)))
                roiIm, roiAff = image_utils.MRIread(g[i])
                roiDef, _ = image_utils.deform(roiIm, roiAff, RASres, device, mode='linear')
                roiDef *= REFmaskRes
                image_utils.MRIwrite(roiDef, affnew, output_directory + '/nonlinearly_registered_roi_' + os.path.split(g[i])[1])

    if Pmesh is not None:
        cmd += (' -f ' + output_directory + '/registered.surf')
        PmeshRot = (np.concatenate([Pmesh[:nv_orig,:], np.ones([nv_orig, 1])], axis=1) @ Tmesh.T)[:,:-1]
        write_geometry(output_directory + '/registered.surf', PmeshRot, TRImesh, volume_info=meta_mesh)

    if arguments.deform_recon_dir is not None:
        print('MRI recon directory provided; working on FreeSurfer data')
        image_utils.deform_FS_derivatives(arguments.deform_recon_dir, arguments.hemisphere, RASres, affnew, output_directory, device)
        cmd += (' -v ' + output_directory + '/aparc+aseg.deformed.mgz:colormap=lut')
        if (arguments.hemisphere=='both') or (arguments.hemisphere=='left'):
            cmd += (' -f ' + output_directory + '/lh.pial.deformed  -f ' + output_directory + '/lh.white.deformed')
        if (arguments.hemisphere=='both') or (arguments.hemisphere=='right'):
            cmd += (' -f ' + output_directory + '/rh.pial.deformed  -f ' + output_directory + '/rh.white.deformed')

    now2 = datetime.now()
    current_time = now2.strftime("%H:%M:%S")
    print("Current Time =", current_time)
    runtime = now2 - now
    print("Running Time =", runtime)

    print(cmd)

# execute script
if __name__ == '__main__':
    main()
