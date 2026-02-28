# Imports
import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as nnf
import ext.interpol as interpol
from photo_reconstruction.image_utils import svf_integration

class photo_aligner(nn.Module):
    """
    Main class to perform alignment
    """

    def __init__(
        self,
        photo_vol,
        photo_dist_vol,
        photo_aff,
        mri_vol,
        mri_mask,
        mri_aff,
        Pmesh,
        Dmesh=None,
        pixel_size=1.0,
        t_ini=None,
        theta_ini=None,
        shear_ini=None,
        scaling_ini=None,
        sz_ini=None,
        y_shifts=None,
        field_ini=None,
        t_mri_ini=None,
        theta_mri_ini=None,
        shear_mri_ini=None,
        s_mri_ini=None,
        field3d_ini=None,
        t_mesh_ini=None,
        theta_mesh_ini=None,
        allow_scaling_and_shear=False,
        allow_sz=False,
        allow_nonlin=False,
        allow_affine_mri=False,
        allow_nonlin_mri=False,
        svf_mode_photos=False,
        cp_spacing2d=None,
        cp_spacing3d=None,
        k_lncc_mri=1.0,
        k_dice_mri=1.0,
        k_dif_slice_loss=1.0,
        k_mesh_loss = 1.0,
        k_regularizer=0.0003,
        k_regularizer_nonlin=0.0003,
        k_regularizer_nonlin3d=0.0001,
        k_regularizer_sz=0.001,
        pad_ignore=None,
        device='cpu'
    ):

        super().__init__()

        # Main data variables
        self.device = device
        self.photo_vol = torch.Tensor(photo_vol.copy()).to(self.device)
        self.photo_rearranged = torch.unsqueeze(torch.unsqueeze(self.photo_vol, dim=0), dim=0).to(self.device)
        self.photo_dist = torch.Tensor(photo_dist_vol.copy()).to(self.device)
        self.photo_dist_rearranged = torch.unsqueeze(torch.unsqueeze(self.photo_dist, dim=0), dim=0).to(self.device)
        self.min_photo_dist = torch.min(self.photo_dist)
        self.photo_aff = torch.Tensor(photo_aff).to(self.device)
        self.svf_mode_photos = svf_mode_photos
        if mri_vol is None:
            self.mri_vol = self.mri_rearranged = self.mri_mask = self.mri_mask_rearranged = self.mri_aff = None
        else:
            self.mri_vol = torch.Tensor(mri_vol).to(self.device)
            self.mri_rearranged = torch.unsqueeze(torch.unsqueeze(self.mri_vol, dim=0), dim=0).to(self.device)
            self.mri_mask = torch.Tensor(mri_mask).to(self.device)
            self.mri_mask_rearranged = torch.unsqueeze(torch.unsqueeze(self.mri_mask, dim=0), dim=0).to(self.device)
            self.mri_aff = torch.Tensor(mri_aff).to(self.device)
        if Pmesh is None:
            self.Pmesh = self.Dmesh = self.Wmesh = None
        else:
            self.Pmesh = torch.Tensor(Pmesh.copy().T).to(self.device)
            self.Dmesh = torch.Tensor(Dmesh.copy()).to(self.device)
            n_v_orig = int(torch.where(self.Dmesh>0)[0][0])
            rel_weight_explorers = n_v_orig / (len(self.Dmesh) - n_v_orig)
            self.Wmesh = torch.concatenate([torch.ones(n_v_orig, device=device), rel_weight_explorers * torch.ones(len(self.Dmesh) - n_v_orig, device=device)])

        # Some constants we'll reuse
        self.photo_siz = self.photo_vol.shape[:-1]
        self.Nslices = self.photo_vol.shape[-1]
        self.pad_ignore = pad_ignore
        self.pixel_size = pixel_size
        self.DELTA = 0.1 / np.log(photo_vol.shape[1])
        self.y_shifts = None if (y_shifts is None) else torch.Tensor(y_shifts).to(self.device)

        # Loss constants
        self.k_lncc_mri = k_lncc_mri
        self.k_dice_mri = k_dice_mri
        self.k_dif_slice_loss = k_dif_slice_loss
        self.k_mesh_loss = k_mesh_loss
        self.k_regularizer = k_regularizer
        self.k_regularizer_nonlin = k_regularizer_nonlin
        self.k_regularizer_nonlin3d = k_regularizer_nonlin3d
        self.k_regularizer_sz = k_regularizer_sz

        # Photo parameters
        if t_ini is not None:
            self.t = torch.nn.Parameter(torch.tensor(t_ini).to(self.device))
        else:
            self.t = torch.nn.Parameter(torch.zeros(2, self.Nslices).to(self.device))
        self.t.requires_grad = True

        if theta_ini is not None:
            self.theta = torch.nn.Parameter(torch.tensor(theta_ini).to(self.device))
        else:
            self.theta = torch.nn.Parameter(torch.zeros(self.Nslices).to(self.device))
        self.theta.requires_grad = True

        if allow_scaling_and_shear:
            if shear_ini is not None:
                self.shear = torch.nn.Parameter(torch.tensor(shear_ini).to(self.device))
            else:
                self.shear = torch.nn.Parameter(torch.zeros(2, self.Nslices).to(self.device))
            self.shear.requires_grad = True
            if scaling_ini is not None:
                self.scaling = torch.nn.Parameter(torch.tensor(scaling_ini).to(self.device))
            else:
                self.scaling = torch.nn.Parameter(torch.zeros(2, self.Nslices).to(self.device))
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
                self.sz = torch.nn.Parameter(torch.tensor(sz_ini).to(self.device))
            else:
                self.sz = torch.nn.Parameter(torch.zeros(1).to(self.device))
            self.sz.requires_grad = True
        else:
            if sz_ini is not None:
                self.sz = torch.tensor(sz_ini).to(self.device)
            else:
                self.sz = torch.zeros(1).to(self.device)

        if allow_nonlin:
            if field_ini is not None:
                self.field = torch.nn.Parameter(torch.tensor(field_ini).to(self.device))
            else:
                ncp = np.ceil(np.array(self.photo_siz) * pixel_size / cp_spacing2d).astype(int)
                self.field = torch.nn.Parameter(torch.zeros(2, ncp[0], ncp[1], self.Nslices).to(self.device))
            self.field.requires_grad = True
        else:
            if field_ini is not None:
                self.field = torch.tensor(field_ini).to(self.device)
            else:
                self.field = None

        # MRI parameters
        if self.mri_vol is None:
            self.t_mri = self.theta_mri = self.shear_mri = self.s_mri = self.field3d = None
        else:
            if t_mri_ini is not None:
                self.t_mri = torch.nn.Parameter(torch.tensor(t_mri_ini).to(self.device))
            else:
                self.t_mri = torch.nn.Parameter(torch.zeros(3).to(self.device))
            self.t_mri.requires_grad = True

            if theta_mri_ini is not None:
                self.theta_mri = torch.nn.Parameter(torch.tensor(theta_mri_ini).to(self.device))
            else:
                self.theta_mri = torch.nn.Parameter(torch.zeros(3).to(self.device))
            self.theta_mri.requires_grad = True

            if allow_affine_mri:
                if shear_mri_ini is not None:
                    self.shear_mri = torch.nn.Parameter(torch.tensor(shear_mri_ini).to(self.device))
                else:
                    self.shear_mri = torch.nn.Parameter(torch.zeros(3).to(self.device))
                self.shear_mri.requires_grad = True
                if s_mri_ini is not None:
                    self.s_mri = torch.nn.Parameter(torch.tensor(s_mri_ini).to(self.device))
                else:
                    self.s_mri = torch.nn.Parameter(torch.zeros(3).to(self.device))
                self.s_mri.requires_grad = True
            else:
                if shear_mri_ini is not None:
                    self.shear_mri = torch.tensor(shear_mri_ini).to(self.device)
                else:
                    self.shear_mri = torch.zeros(3).to(self.device)
                if s_mri_ini is not None:
                    self.s_mri = torch.tensor(s_mri_ini).to(self.device)
                else:
                    self.s_mri = torch.zeros(3).to(self.device)

            if allow_nonlin_mri:
                if field3d_ini is not None:
                    self.field3d = torch.nn.Parameter(torch.tensor(field3d_ini).to(self.device))
                else:
                    ncp = np.ceil(np.array(self.mri_vol.shape) / cp_spacing3d).astype(int)
                    self.field3d = torch.nn.Parameter(torch.zeros(3, ncp[0], ncp[1], ncp[2]).to(self.device))
                self.field3d.requires_grad = True
            else:
                if field3d_ini is not None:
                    self.field3d = torch.tensor(field3d_ini).to(self.device)
                else:
                    self.field3d = None

        # Mesh parameters
        if self.Pmesh is None:
            self.t_mesh = self.theta_mesh = None
        else:
            if t_mesh_ini is not None:
                self.t_mesh = torch.nn.Parameter(torch.tensor(t_mesh_ini).to(self.device))
            else:
                self.t_mesh = torch.nn.Parameter(torch.zeros(3).to(self.device))
            self.t_mesh.requires_grad = True
            if theta_mesh_ini is not None:
                self.theta_mesh = torch.nn.Parameter(torch.tensor(theta_mesh_ini).to(self.device))
            else:
                self.theta_mesh = torch.nn.Parameter(torch.zeros(3).to(self.device))
            self.theta_mesh.requires_grad = True

        # create sampling grid we'll reuse over and over
        vectors = [torch.arange(0, s) for s in self.photo_vol.shape]
        self.grids = torch.stack(torch.meshgrid(vectors)).to(self.device)


    ##########################################
    # Functions to get subsets of parameters #
    ##########################################

    def parameters_photos(self):
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
    def parameters_mri(self):
        if isinstance(self.t_mri, nn.Parameter):
            yield self.t_mri
        if isinstance(self.theta_mri, nn.Parameter):
            yield self.theta_mri
        if isinstance(self.shear_mri, nn.Parameter):
            yield self.shear_mri
        if isinstance(self.s_mri, nn.Parameter):
            yield self.s_mri
        if isinstance(self.field3d, nn.Parameter):
            yield self.field3d
    def parameters_mesh(self):
        if isinstance(self.t_mesh, nn.Parameter):
            yield self.t_mesh
        if isinstance(self.theta_mesh, nn.Parameter):
            yield self.theta_mesh


    ############################
    # The big forward function #
    ############################

    def forward(self):
        # We scale angles / shearings / scalings as a simple form of preconditioning (which shouldn't be needed with bfgs, but whatever...)

        # Parameters of photos
        theta_f = ((self.theta - self.theta.mean()) / 180 * torch.tensor(np.pi)) # subtracting mean kills solution ambiguity with rotation
        shear_f = self.shear / 100  # percentages
        scaling_f = torch.exp(self.scaling / 20)  # ensures positive and symmetry around 1 in log scale
        t_f = self.t + self.DELTA  # no scaling
        sz_f = torch.exp(self.sz / 20)  # ensures positive and symmetry around 1 in log scale
        if self.field is not None:
            # make it a % of the dimensions
            field_x = self.field[0, :, :, :] * torch.Tensor(np.array(self.photo_siz[0] / 100.0) )
            field_y = self.field[1, :, :, :] * torch.Tensor(np.array(self.photo_siz[1] / 100.0) )
            field_pixels = torch.stack([field_x, field_y])
        else:
            field_pixels = None

        # Parameters of reference volume
        if self.mri_vol is None:
            theta_mri_f = t_mri_f = s_mri_f = shear_mri_f = field3d_pixels = None
        else:
            theta_mri_f = (self.theta_mri / 180 * torch.tensor(np.pi))  # degrees -> radians
            t_mri_f = self.t_mri  # no scaling
            s_mri_f = torch.exp(self.s_mri / 20)  # ensures positive and symmetry around 1 in log scale
            shear_mri_f = self.shear_mri / 100 # percentages
            if self.field3d is not None:
                # make it a % of the dimensions
                field3d_x = self.field3d[0, :, :, :] * torch.Tensor(np.array(self.mri_vol.shape[0] / 100.0))
                field3d_y = self.field3d[1, :, :, :] * torch.Tensor(np.array(self.mri_vol.shape[1] / 100.0))
                field3d_z = self.field3d[2, :, :, :] * torch.Tensor(np.array(self.mri_vol.shape[2] / 100.0))
                field3d_pixels = torch.stack([field3d_x, field3d_y, field3d_z])
            else:
                field3d_pixels = None # no non-linear deformation of the MRI

        # Parameters of reference mesh
        if self.Pmesh is None:
            theta_mesh_f = t_mesh_f = None
        else:
            theta_mesh_f = (self.theta_mesh / 180 * torch.tensor(np.pi))  # degrees -> radians
            t_mesh_f = self.t_mesh

        # Prepare 2D matrices for the photos
        M = torch.zeros([4, 4, self.Nslices]).to(self.device)
        M[0, 0, :] = scaling_f[0, :] * (torch.cos(theta_f) - shear_f[0, :] * torch.sin(theta_f))
        M[0, 2, :] = scaling_f[0, :] * (shear_f[1, :] * torch.cos(theta_f) - (1 + shear_f[0, :] * shear_f[1, :]) * torch.sin(theta_f))
        M[0, 3, :] = t_f[0, :]
        M[1, 1, :] = 1
        M[2, 0, :] = scaling_f[1, :] * (torch.sin(theta_f) + shear_f[0, :] * torch.cos(theta_f))
        M[2, 2, :] = scaling_f[1, :] * (shear_f[1, :] * torch.sin(theta_f) + (1 + shear_f[0, :] * shear_f[1, :]) * torch.cos(theta_f))
        M[2, 3, :] = t_f[1, :]
        M[3, 3, :] = 1

        # update mesh grids for photos
        # First, upscale field
        if self.field is not None:
            # Bspline better than trilinear
            # field_fullsiz = torch.nn.Upsample(size=self.photo_siz,  align_corners=True, mode="bilinear" )(field_pixels.permute([0,3,1,2]))
            field_fullsiz = interpol.resize(field_pixels.permute([0,3,1,2]), shape=self.photo_siz, anchor = 'e', interpolation = 3, bound = 'dft', prefilter = False)
            field_fullsiz_rearranged = field_fullsiz.permute(0, 2, 3, 1)
            if self.svf_mode_photos: # integrate field if needed
                field_fullsiz_rearranged = svf_integration(field_fullsiz_rearranged)
            # Membrane energy
            grad_x = (field_fullsiz_rearranged[:, 2:, 1:-1, :] - field_fullsiz_rearranged[:, :-2, 1:-1, :]) / (2.0 * self.pixel_size)
            grad_y = (field_fullsiz_rearranged[:, 1:-1, 2:, :] - field_fullsiz_rearranged[:, 1:-1, :-2, :]) / (2.0 * self.pixel_size)
            cost_field = torch.sum(grad_x * grad_x + grad_y * grad_y) / grad_x[0].numel()
        else:
            field_fullsiz = None
            field_fullsiz_rearranged = None
            cost_field = torch.zeros(1).to(self.device)

        # update mesh grids for photos
        photo_aff = torch.zeros(4, 4).to(self.device)
        photo_aff[:, :] = self.photo_aff
        photo_aff[1, 2:] = self.photo_aff[1, 2:] * sz_f
        if field_fullsiz_rearranged is None:
            nonlingrid = self.grids
        else:
            nonlingrid = self.grids + torch.cat([field_fullsiz_rearranged, torch.zeros_like(self.grids[0:1])], dim=0)
        T = []
        for z in range(self.Nslices):
            T.append(torch.matmul(torch.matmul(torch.inverse(photo_aff), M[:, :, z]), photo_aff ))
        T = torch.stack(T, dim=-1)
        grids_new = []
        for d in range(3):
            grids_new.append(T[d, 0, :][None, None] * nonlingrid[0]+ T[d, 1, :][None, None] * nonlingrid[1] \
                             + T[d, 2, :][None, None] * nonlingrid[2] + T[d, 3, :][None, None])
        grids_new = torch.stack(grids_new, dim=0)

        # Resample photos
        # We need to make the new grid compatible with grid_resample...
        grids_new = torch.unsqueeze(grids_new, 0)
        grids_new = grids_new.permute(0, 2, 3, 4, 1)
        for i in range(3):
            grids_new[:, :, :, :, i] = 2 * (
                grids_new[:, :, :, :, i] / (self.photo_vol.shape[i] - 1)
                - 0.5
            )
        # Not sure why, but channels need to be reversed
        grids_new = grids_new[..., [2, 1, 0]]
        photo_resampled = nnf.grid_sample(self.photo_rearranged, grids_new, align_corners=True, mode="bilinear", padding_mode="zeros")
        photo_resampled = torch.squeeze(photo_resampled.permute(2, 3, 4, 1, 0))
        # Careful: border extrapolation for the disntace maps!
        photo_dist_resampled = nnf.grid_sample(self.photo_dist_rearranged, grids_new, align_corners=True, mode='bilinear', padding_mode='border')
        photo_dist_resampled = torch.squeeze(photo_dist_resampled.permute(2, 3, 4, 1, 0))

        # Now work on the reference
        cost_field3d = torch.Tensor([0.0]).to(self.device)
        mri_aff_combined = Rt = grids_new_mri_nonlin = None
        if self.mri_vol is not None:
            Rx = torch.zeros([4, 4]).to(self.device)
            Rx[0, 0] = 1
            Rx[1, 1] = torch.cos(theta_mri_f[0])
            Rx[1, 2] = -torch.sin(theta_mri_f[0])
            Rx[2, 1] = torch.sin(theta_mri_f[0])
            Rx[2, 2] = torch.cos(theta_mri_f[0])
            Rx[3, 3] = 1

            Ry = torch.zeros([4, 4]).to(self.device)
            Ry[0, 0] = torch.cos(theta_mri_f[1])
            Ry[0, 2] = torch.sin(theta_mri_f[1])
            Ry[1, 1] = 1
            Ry[2, 0] = -torch.sin(theta_mri_f[1])
            Ry[2, 2] = torch.cos(theta_mri_f[1])
            Ry[3, 3] = 1

            Rz = torch.zeros([4, 4]).to(self.device)
            Rz[0, 0] = torch.cos(theta_mri_f[2])
            Rz[0, 1] = -torch.sin(theta_mri_f[2])
            Rz[1, 0] = torch.sin(theta_mri_f[2])
            Rz[1, 1] = torch.cos(theta_mri_f[2])
            Rz[2, 2] = 1
            Rz[3, 3] = 1

            SHx = torch.zeros([4, 4]).to(self.device)
            SHx[0, 0] = 1
            SHx[1, 0] = shear_mri_f[1]
            SHx[1, 1] = 1
            SHx[2, 0] = shear_mri_f[2]
            SHx[2, 2] = 1
            SHx[3, 3] = 1

            SHy = torch.zeros([4, 4]).to(self.device)
            SHy[0, 0] = 1
            SHy[0, 1] = shear_mri_f[0]
            SHy[1, 1] = 1
            SHy[2, 1] = shear_mri_f[2]
            SHy[2, 2] = 1
            SHy[3, 3] = 1

            # SHz = np.array([[1, 0, sh[0]], [0, 1, sh[1]], [0, 0, 1]])
            SHz = torch.zeros([4, 4]).to(self.device)
            SHz[0, 0] = 1
            SHz[0, 2] = shear_mri_f[0]
            SHz[1, 1] = 1
            SHz[1, 2] = shear_mri_f[1]
            SHz[2, 2] = 1
            SHz[3, 3] = 1

            trans_and_scale = torch.eye(4).to(self.device)
            trans_and_scale[:-1, -1] = t_mri_f
            trans_and_scale[0, 0] = s_mri_f[0]
            trans_and_scale[1, 1] = s_mri_f[1]
            trans_and_scale[2, 2] = s_mri_f[2]

            Rt = trans_and_scale @ SHx @ SHy @ SHz @ Rx @ Ry @ Rz

            # Now let's work on the deformed volume
            mri_aff_combined = torch.matmul(Rt, self.mri_aff)
            if self.y_shifts is None:
                # without shifts, it's all linear in one step
                grids_new_mri = torch.zeros(self.grids.shape).to(self.device)
                D = torch.matmul(torch.inverse(mri_aff_combined), photo_aff)
                for d in range(3):
                    grids_new_mri[d, :, :, :] = (D[d, 0] * self.grids[0, :, :, :]
                        + D[d, 1] * self.grids[1, :, :, :]
                        + D[d, 2] * self.grids[2, :, :, :]
                        + grids_new_mri[d, :, :, :]
                        + D[d, 3]
                    )
            else:
                # with shifts: we do RAS first, shift A-P, and then go to vox in MRI
                grids_photo_ras = torch.zeros(self.grids.shape).to(self.device)
                for d in range(3):
                    grids_photo_ras[d, :, :, :] = (photo_aff[d, 0] * self.grids[0, :, :, :]
                        + photo_aff[d, 1] * self.grids[1, :, :, :]
                        + photo_aff[d, 2] * self.grids[2, :, :, :]
                        + grids_photo_ras[d, :, :, :]
                        + photo_aff[d, 3]
                    )
                grids_photo_ras[1] -= (sz_f * self.y_shifts[None, None, :])  # crucial to scale shifts
                grids_new_mri = torch.zeros(self.grids.shape).to(self.device)
                D = torch.inverse(mri_aff_combined)
                for d in range(3):
                    grids_new_mri[d, :, :, :] = (D[d, 0] * grids_photo_ras[0, :, :, :]
                        + D[d, 1] * grids_photo_ras[1, :, :, :]
                        + D[d, 2] * grids_photo_ras[2, :, :, :]
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

            # Resample field if needed
            if self.field3d is  None:
                field3d_fullsiz_rearranged = None
                grids_new_mri_nonlin = grids_new_mri
            else:
                # Upscale field and compute membrane energy TODO: no need to go to full size, can interpolate directly in low-res (see PSAMSEG code)
                # Replace trilinear by Bspline
                # field3d_fullsiz = torch.nn.Upsample(size=self.mri_vol.shape, align_corners=True, mode="trilinear")(torch.unsqueeze(field3d_pixels, dim=0))
                field3d_fullsiz = interpol.resize(torch.unsqueeze(field3d_pixels, dim=0), shape=self.mri_vol.shape, anchor = 'e', interpolation = 3, bound = 'dft',  prefilter = False)
                field3d_fullsiz_rearranged = torch.squeeze(field3d_fullsiz)
                x = (field3d_fullsiz_rearranged[:, 2:, 1:-1, 1:-1] - field3d_fullsiz_rearranged[:, :-2, 1:-1, 1:-1]) / (2.0)
                y = (field3d_fullsiz_rearranged[:, 1:-1, 2:, 1:-1] - field3d_fullsiz_rearranged[:, 1:-1, :-2, 1:-1]) / (2.0)
                z = (field3d_fullsiz_rearranged[:, 1:-1, 1:-1, 2:] - field3d_fullsiz_rearranged[:, 1:-1, 1:-1, :-2]) / (2.0)
                cost_field3d = torch.sum(x*x + y*y + z*z) / field3d_fullsiz_rearranged[0].numel()

                # Interpolate shifts and add to deformation
                field_resampled = nnf.grid_sample(field3d_fullsiz, grids_new_mri, align_corners=True, mode="bilinear", padding_mode="zeros")
                field_resampled_rearranged = field_resampled.permute([0,2,3,4,1])[..., [2, 1, 0]]
                grids_new_mri_nonlin = torch.zeros_like(grids_new_mri, device=self.device)
                for d in range(3):
                    grids_new_mri_nonlin[:, :, :, :, d] = grids_new_mri[:, :, :, :, d] + 2 * (field_resampled_rearranged[:, :, :, :, d] / (self.mri_vol.shape[2-d] - 1))

            # acual deformation of  MRI and corresponding mask
            mri_resampled = nnf.grid_sample(self.mri_rearranged, grids_new_mri_nonlin, align_corners=True, mode="bilinear", padding_mode="zeros")
            mri_resampled = torch.squeeze(mri_resampled)
            mri_mask_resampled = nnf.grid_sample(self.mri_mask_rearranged, grids_new_mri_nonlin, align_corners=True, mode="bilinear", padding_mode="zeros")
            mri_mask_resampled = torch.squeeze(mri_mask_resampled)

        # Let's work on the mesh coordinates, if available
        if self.Pmesh is None:
            Tmesh = None
        else:
            RxM = torch.zeros([4, 4]).to(self.device)
            RxM[0, 0] = 1
            RxM[1, 1] = torch.cos(theta_mesh_f[0])
            RxM[1, 2] = -torch.sin(theta_mesh_f[0])
            RxM[2, 1] = torch.sin(theta_mesh_f[0])
            RxM[2, 2] = torch.cos(theta_mesh_f[0])
            RxM[3, 3] = 1
            RyM = torch.zeros([4, 4]).to(self.device)
            RyM[0, 0] = torch.cos(theta_mesh_f[1])
            RyM[0, 2] = torch.sin(theta_mesh_f[1])
            RyM[1, 1] = 1
            RyM[2, 0] = -torch.sin(theta_mesh_f[1])
            RyM[2, 2] = torch.cos(theta_mesh_f[1])
            RyM[3, 3] = 1
            RzM = torch.zeros([4, 4]).to(self.device)
            RzM[0, 0] = torch.cos(theta_mesh_f[2])
            RzM[0, 1] = -torch.sin(theta_mesh_f[2])
            RzM[1, 0] = torch.sin(theta_mesh_f[2])
            RzM[1, 1] = torch.cos(theta_mesh_f[2])
            RzM[2, 2] = 1
            RzM[3, 3] = 1
            trans_and_scaleM = torch.eye(4).to(self.device)
            trans_and_scaleM[:-1, -1] = t_mesh_f
            Tmesh = trans_and_scaleM @ RxM @ RyM @ RzM
            # Computing coordinates in vox space
            if self.y_shifts is None: # simple matrix multiplications if there's no shifts...
                RtM = torch.inverse(photo_aff) @  Tmesh
                vox = RtM @ torch.concat([self.Pmesh, torch.ones([1, self.Pmesh.shape[1]],device=self.device)], dim=0)
            else: # otherwise we need to be careful to check in which "sandwich" we found ourselves!
                PmeshTrans = Tmesh @ torch.concat([self.Pmesh, torch.ones([1, self.Pmesh.shape[1]], device=self.device)], dim=0)
                photo_aff_inv = torch.inverse(photo_aff)
                i_coords = photo_aff_inv[0, 0] * PmeshTrans[0] + photo_aff_inv[0, 2] * PmeshTrans[2] + photo_aff_inv[0, 3]
                j_coords = photo_aff_inv[1, 0] * PmeshTrans[0] + photo_aff_inv[1, 2] * PmeshTrans[2] + photo_aff_inv[1, 3]
                o_coords = torch.ones_like(i_coords)
                yp = photo_aff[1, 3] + photo_aff[1, 2] * torch.arange(self.Nslices, device=self.device) - self.y_shifts * sz_f # crucial to scale shifts
                k_coords = -torch.ones_like(i_coords)
                for z in range(self.Nslices-1):
                    ok = (PmeshTrans[1]<=yp[z]) & (PmeshTrans[1]>yp[z+1])
                    da = (yp[z] - PmeshTrans[1][ok]).clip(1e-6)
                    db = (PmeshTrans[1][ok] - yp[z+1]).clip(1e-6)
                    wa = db / (da + db)
                    wb = 1 - wa
                    k_coords[ok] = wa * z + wb * (z+1)
                vox = torch.stack([i_coords, j_coords, k_coords, o_coords])

            ok = ( (vox[0,:] > 0) & (vox[1,:] > 0) & (vox[2,:] > 0) & (vox[0,:] <= photo_dist_resampled.shape[0] - 1)
                       & (vox[1,:] <= photo_dist_resampled.shape[1] - 1) & (vox[2,:] <= photo_dist_resampled.shape[2] - 1) )
            IIv = vox[0, ok]
            JJv = vox[1, ok]
            KKv = vox[2, ok]
            fx = torch.floor(IIv).long()
            cx = fx + 1
            cx[cx > (photo_dist_resampled.shape[0] - 1)] = (photo_dist_resampled.shape[0] - 1)
            wcx = IIv - fx
            wfx = 1 - wcx
            fy = torch.floor(JJv).long()
            cy = fy + 1
            cy[cy > (photo_dist_resampled.shape[1] - 1)] = (photo_dist_resampled.shape[1] - 1)
            wcy = JJv - fy
            wfy = 1 - wcy
            fz = torch.floor(KKv).long()
            cz = fz + 1
            cz[cz > (photo_dist_resampled.shape[2] - 1)] = (photo_dist_resampled.shape[2] - 1)
            wcz = KKv - fz
            wfz = 1 - wcz
            c000 = photo_dist_resampled[fx, fy, fz]
            c100 = photo_dist_resampled[cx, fy, fz]
            c010 = photo_dist_resampled[fx, cy, fz]
            c110 = photo_dist_resampled[cx, cy, fz]
            c001 = photo_dist_resampled[fx, fy, cz]
            c101 = photo_dist_resampled[cx, fy, cz]
            c011 = photo_dist_resampled[fx, cy, cz]
            c111 = photo_dist_resampled[cx, cy, cz]
            c00 = c000 * wfx + c100 * wcx
            c01 = c001 * wfx + c101 * wcx
            c10 = c010 * wfx + c110 * wcx
            c11 = c011 * wfx + c111 * wcx
            c0 = c00 * wfy + c10 * wcy
            c1 = c01 * wfy + c11 * wcy
            values = c0 * wfz + c1 * wcz
            # New symmetric cost
            # cost_mesh = torch.sum(torch.abs(values) * self.Wmesh[ok]) / torch.sum(self.Wmesh[ok])
            x = torch.abs(values)
            y = self.Dmesh[ok]
            rho = 1.0
            w_x = torch.exp(-rho * x)
            w_y = torch.exp(-rho * y) * self.Wmesh[ok]
            cost_mesh = torch.sum(w_x * y) / w_x.sum() + torch.sum(w_y * x) / w_y.sum()


        # now let's compute the metrics
        # First, LNCC photos<->MRI (computed in 2D, of course, slice thickness is too large)
        if self.mri_vol is None:
            lncc_mri_loss = torch.zeros(1).to(self.device)
            mri_resampled_masked = None
        else:
            mri_resampled_masked = mri_resampled * mri_mask_resampled
            Ii = mri_resampled_masked[:, :, self.pad_ignore: -self.pad_ignore].permute([2, 0, 1])[:, None, :, :]
            Ji = photo_resampled[:, :, self.pad_ignore : -self.pad_ignore].permute([2,0,1])[:,None,:,:] / 255.0
            sum_filt = torch.ones([1, 1, 9, 9]).to(self.device)
            stride = (1, 1)
            padding = (4,4)
            I2 = Ii * Ii; J2 = Ji * Ji; IJ = Ii * Ji
            I_sum = torch.nn.functional.conv2d(Ii, sum_filt, stride=stride, padding=padding)
            J_sum = torch.nn.functional.conv2d(Ji, sum_filt, stride=stride, padding=padding)
            I2_sum = torch.nn.functional.conv2d(I2, sum_filt, stride=stride, padding=padding)
            J2_sum = torch.nn.functional.conv2d(J2, sum_filt, stride=stride, padding=padding)
            IJ_sum = torch.nn.functional.conv2d(IJ, sum_filt, stride=stride, padding=padding)
            u_I = I_sum / 81.0
            u_J = J_sum / 81.0
            cross = IJ_sum - u_J * I_sum - u_I * J_sum + u_I * u_J * 81.0
            I_var = I2_sum - 2 * u_I * I_sum + u_I * u_I * 81.0
            J_var = J2_sum - 2 * u_J * J_sum + u_J * u_J * 81.0
            cc = cross * cross / (I_var * J_var + 1e-5)
            lncc_mri_loss = 1 - torch.mean(cc)

        # Next: Dice between photo masks and MRI mask
        if self.mri_vol is None:
            dice_mri_loss = torch.zeros(1).to(self.device)
        else:
            rho = 2.0
            A = mri_mask_resampled
            B = torch.sigmoid(rho * photo_dist_resampled)
            dice = 2.0 * torch.mean(A * B) / (1e-6 + torch.mean(A * A) + torch.mean(B * B))
            dice_mri_loss = 1 - dice

        # Next: mean absolute difference between consecutive slices
        dif = (photo_resampled[:, :, self.pad_ignore + 1: - self.pad_ignore] -
            photo_resampled[:, :, self.pad_ignore: - self.pad_ignore - 1])
        if True:  # maybe give double weight to first and last? The other slices count twice after all
            dif[:, :, 1] *= 2.0
            dif[:, :, -1] *= 2.0
        dif_slice_loss = torch.mean(torch.abs(dif)) / 255

        # The last data term is the surface
        if self.Pmesh is None:
            mesh_loss = torch.zeros(1).to(self.device)
        else:
            mesh_loss = cost_mesh * 0.1 # we scale more lor less to 0-1

        # Finally, the affine regularizer
        aff_regularizers = torch.abs(torch.sum(self.scaling, dim=0) / 20)
        loss_aff_regularizer = torch.mean(aff_regularizers)

        if False:
            from photo_reconstruction.image_utils import MRIread, MRIwrite
            MRIwrite(mri_resampled_masked.cpu().detach().numpy(), photo_aff.cpu().detach().numpy(), '/tmp/kk1.mgz')
            MRIwrite(photo_resampled.cpu().detach().numpy(), photo_aff.cpu().detach().numpy(), '/tmp/kk2.mgz')
            bif = torch.zeros_like(mri_resampled_masked)
            bif[:,:, self.pad_ignore : - self.pad_ignore] = cc[:,0].permute([1,2,0])
            MRIwrite(bif.cpu().detach().numpy(), photo_aff.cpu().detach().numpy(), '/tmp/kk3.mgz')
            bif = torch.zeros_like(mri_resampled_masked)
            vox2 = vox.clone(); vox2[vox2<0]=0;
            for kk in range(3):
                vox2[kk,vox2[kk,:]>(bif.shape[kk]-1)]=(bif.shape[kk]-1)
            vox2 = vox2[:-1,:].long()
            bif[vox2[0,:],vox2[1,:],vox2[2,:]]=1
            MRIwrite(bif.cpu().detach().numpy(), photo_aff.cpu().detach().numpy(), '/tmp/kk4.mgz')
            print('freeview /tmp/kk1.mgz /tmp/kk2.mgz /tmp/kk3.mgz /tmp/kk4.mgz')

        TINY = 1e-6
        loss = (
                self.k_lncc_mri * lncc_mri_loss
                + self.k_dice_mri * dice_mri_loss
                + self.k_dif_slice_loss * dif_slice_loss
                + self.k_mesh_loss * mesh_loss
                + self.k_regularizer * loss_aff_regularizer
                + self.k_regularizer_nonlin * cost_field
                + self.k_regularizer_nonlin3d * cost_field3d
                + self.k_regularizer_sz * self.sz * self.sz
                + TINY * (torch.mean(torch.square((self.t - self.DELTA))))
                + TINY * (torch.mean(torch.square((self.theta))))
                + TINY * (torch.mean(torch.square((self.shear))))
                + TINY * (torch.mean(torch.square((self.scaling))))
        )
        if torch.isnan(loss):
            if True:
                import pdb; pdb.set_trace()
            else:
                raise Exception('NaN in loss...')

        return loss, photo_resampled, photo_aff, mri_aff_combined, Rt, T, mri_resampled_masked, Tmesh, grids_new_mri_nonlin

