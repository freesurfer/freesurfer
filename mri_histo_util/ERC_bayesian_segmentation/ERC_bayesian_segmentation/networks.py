import torch
import torch.nn as nn
import numpy as np
import math
import ext.interpol as interpol
from ext.my_functions import membrane3


class AtlasDeformer(nn.Module):
    """
    Main class to perform alignment
    """

    def __init__(self,
                 I, aff_I, A, aff_A, ts_ini=None, thetas_ini=None, scalings_ini=None,
                 shears_ini=None, mus=None, vars=None, mu_bg=0, var_bg=100, weights=None,
                 gmm_components=None, k_penalty_grad=0.01, FIELD_ini=None, cp_spacing=20, device='cpu',
                 allow_shearing=False, allow_nonlinear=False, allow_similarity=True,
                 skip=1, dtype=torch.float32):

        super().__init__()

        backend = dict(dtype=dtype, device=device)
        self.device = device
        self.dtype = dtype

        self.skip = skip
        self.mu_bg = torch.tensor(np.array(mu_bg), **backend)
        self.var_bg = torch.tensor(np.array(var_bg), **backend)
        self.mus = torch.tensor(mus, **backend)
        self.vars = torch.tensor(vars, **backend)
        self.weights = torch.tensor(weights, **backend)
        self.gmm_components = torch.tensor(gmm_components.astype('int'), **backend)
        self.num_gmm_components = int(self.gmm_components.sum().item())

        self.siz = np.asarray(I.shape, dtype='int')
        self.Iskip = torch.tensor(np.array(I[::skip, ::skip, ::skip]), **backend)
        self.Iskip_siz = self.Iskip.shape
        self.Iskip_siz_tensor = torch.tensor(self.Iskip.shape, **backend)
        self.aff_I = torch.tensor(aff_I, **backend)
        self.A = torch.tensor(A, **backend)
        self.A_rearranged = torch.unsqueeze(self.A.permute(3, 0, 1, 2), dim=0)
        self.sizA = np.asarray(A.shape[:-1], dtype='int')
        self.sizA_tensor = torch.tensor(self.sizA, **backend)
        self.aff_A = torch.tensor(aff_A, **backend)
        self.atlas_voxsize = torch.abs(torch.det(self.aff_A)) ** (1./3.)
        self.cp_spacing = cp_spacing
        self.k_penalty_grad = torch.tensor(np.array(k_penalty_grad), **backend)
        self.n_classes = A.shape[-1]
        self.siz_tensor = torch.tensor(self.siz, **backend)
        self.volres = np.sqrt(np.sum(aff_I[:-1, :-1] ** 2, axis=0))
        self.eps = torch.tensor(np.array(1e-9), **backend)
        self.GAUSSIAN_LHOODS = torch.zeros([*self.Iskip_siz, self.num_gmm_components+1], **backend)
        self.GAUSSIAN_LHOODS[:, :, :, 0] = 1 / torch.sqrt(2 * math.pi * self.var_bg) * torch.exp(
            -0.5 * torch.pow(self.Iskip - self.mu_bg, 2.0) / self.var_bg)
        for c in range(self.num_gmm_components):
            self.GAUSSIAN_LHOODS[:, :, :, c+1] = 1.0 / torch.sqrt(2 * math.pi * self.vars[c]) * torch.exp(
                -0.5 * torch.pow(self.Iskip - self.mus[c], 2.0) / self.vars[c])

        if ts_ini is not None:
            self.ts = torch.tensor(ts_ini, **backend)
        else:
            self.ts = torch.zeros(3, **backend)
        if thetas_ini is not None:
            self.thetas = torch.tensor(thetas_ini, **backend)
        else:
            self.thetas = torch.zeros(3, **backend)
        if scalings_ini is not None:
            self.scalings = torch.tensor(scalings_ini, **backend)
        else:
            self.scalings = torch.zeros(3, **backend)
        if allow_similarity:
            self.ts = torch.nn.Parameter(self.ts, requires_grad=True)
            self.thetas = torch.nn.Parameter(self.thetas, requires_grad=True)
            self.scalings = torch.nn.Parameter(self.scalings, requires_grad=True)

        if shears_ini is not None:
            self.shears = torch.tensor(shears_ini, **backend)
        else:
            self.shears = torch.zeros(3, **backend)
        if allow_shearing:
            self.shears = torch.nn.Parameter(self.shears, requires_grad=True)

        self.FIELD_siz = np.round(self.siz * self.volres / cp_spacing).astype('int')
        if FIELD_ini is not None:
            self.FIELD = torch.zeros([*self.FIELD_siz, 3], **backend)
            factor = self.FIELD_siz / FIELD_ini.shape[:-1]
            FIELD_ini = torch.tensor(FIELD_ini, **backend)
            for c in range(3):
                self.FIELD[:, :, :, c] = interpol.resize(FIELD_ini[:, :, :, c], shape=self.FIELD_siz.tolist(), anchor='e', interpolation=3, bound='dft', prefilter=False)
                self.FIELD[:, :, :, c] = interpol.spline_coeff_nd(self.FIELD[:, :, :, c], bound='dft', interpolation=3)
        else:
            self.FIELD = torch.zeros([*self.FIELD_siz, 3], **backend)
        if allow_nonlinear:
            self.FIELD = torch.nn.Parameter(self.FIELD, requires_grad=True)

        self.scale_term = self.siz_tensor/torch.tensor(self.FIELD_siz, device=device)
        # create sampling grid
        vectors = [torch.arange(0, s, skip, device=device) for s in self.siz]
        self.grids = torch.stack(torch.meshgrid(vectors))

    def forward(self):

        grids_new_atlas, field_f = self.get_atlas_sampling_grid()

        # Theoretically you should be able to scale the penalty by 1/(control_point_reso/image_reso)^2
        # The scale term: (image_size/field_size) is basically just that. Could change to the resolution
        # factor above if we want to. Note btw that the analytical loss expects isotropic resolution.
        # Also image_reso and atlas_reso are the same so doesn't matter which one you use, but it should
        # be the atlas_reso so 1/(control_point_reso/atlas_reso)^2
        grad_penalty = self.get_grad_penalty_new(field_f/self.scale_term.view([1, 1, 1, 3]))
        # Prepare the new grid for interpolation

        atlas_resampled = interpol.grid_pull(self.A_rearranged, grids_new_atlas, interpolation=1)
        atlas_resampled = torch.clamp(torch.squeeze(atlas_resampled.permute(2, 3, 4, 1, 0)), min=0, max=1)

        # # Extrapolation: put prior on background
        missing_mass = 1 - torch.sum(atlas_resampled,dim=3)
        prior_bg = (atlas_resampled[:, :, :, 0] + missing_mass).clamp(0, 1)
        priors_rest = (atlas_resampled[:, :, :, 1:]).clamp(0, 1)
        priors = torch.cat([torch.unsqueeze(prior_bg, dim=-1), priors_rest], dim=-1)

        # # Everything's now in place! We can compute the log posteriors
        LHOODS = torch.zeros_like(priors)
        LHOODS[:, :, :, 0] = prior_bg * self.GAUSSIAN_LHOODS[:, :, :, 0]
        for c in range(self.n_classes - 1):
            num_gaussians = self.gmm_components[c].int()
            for g in range(num_gaussians):
                gaussian_number = sum(self.gmm_components[:c].int()) + g
                LHOODS[:, :, :, c + 1] += self.weights[gaussian_number] * priors_rest[:, :, :, c] * self.GAUSSIAN_LHOODS[:, :, :, gaussian_number + 1]

        # LHOODS, priors = self.get_priors(grids_new_atlas)
        normalizer = torch.sum(LHOODS,dim=-1)
        LOG_LHOOD = torch.log(self.eps+normalizer)

        deformation_penalty = self.k_penalty_grad * grad_penalty
        loss = - torch.mean(LOG_LHOOD) + deformation_penalty

        if torch.isnan(loss):
            raise Exception('nan in loss...')

        return loss, grids_new_atlas, priors, deformation_penalty


    def get_atlas_sampling_grid(self):

        # We scale angles / shearings / scalings as a simple form of precondit  ioning
        thetas_f = self.thetas / 180 * math.pi  # degrees -> radians
        shears_f = self.shears / 100  # percentages
        scalings_f = torch.exp(self.scalings / 20)  # ensures positive and symmetry around 1 in log scale
        ts_f = self.ts  # in mm
        field_f = self.FIELD / 100 * self.sizA_tensor.view([1,1,1,3])  # percentages

        # Prepare affine matrix for the atlas deformation
        # Rotations matrices
        Rx = torch.zeros([4, 4], dtype=self.dtype, device=self.device)
        Rx[0, 0] = 1
        Rx[1, 1] = torch.cos(thetas_f[0])
        Rx[1, 2] = -torch.sin(thetas_f[0])
        Rx[2, 1] = torch.sin(thetas_f[0])
        Rx[2, 2] = torch.cos(thetas_f[0])
        Rx[3, 3] = 1

        Ry = torch.zeros([4, 4], dtype=self.dtype, device=self.device)
        Ry[0, 0] = torch.cos(thetas_f[1])
        Ry[0, 2] = torch.sin(thetas_f[1])
        Ry[1, 1] = 1
        Ry[2, 0] = -torch.sin(thetas_f[1])
        Ry[2, 2] = torch.cos(thetas_f[1])
        Ry[3, 3] = 1

        Rz = torch.zeros([4, 4], dtype=self.dtype, device=self.device)
        Rz[0, 0] = torch.cos(thetas_f[2])
        Rz[0, 1] = -torch.sin(thetas_f[2])
        Rz[1, 0] = torch.sin(thetas_f[2])
        Rz[1, 1] = torch.cos(thetas_f[2])
        Rz[2, 2] = 1
        Rz[3, 3] = 1

        # Scaling matrix
        Sc = torch.zeros([4, 4], dtype=self.dtype, device=self.device)
        Sc[0, 0] = scalings_f[0]
        Sc[1, 1] = scalings_f[1]
        Sc[2, 2] = scalings_f[2]
        Sc[3, 3] = 1

        # Shearing matrix
        Sh = torch.zeros([4, 4], dtype=self.dtype, device=self.device)
        Sh[0, 0] = 1
        Sh[0, 1] = shears_f[1]
        Sh[0, 2] = shears_f[2]
        Sh[1, 0] = shears_f[0]
        Sh[1, 1] = 1
        Sh[1, 2] = shears_f[2]
        Sh[2, 0] = shears_f[0]
        Sh[2, 1] = shears_f[1]
        Sh[2, 2] = 1
        Sh[3, 3] = 1

        # Translation matrix
        T = torch.eye(4, dtype=self.dtype, device=self.device)
        T[0, 0] = 1
        T[1, 1] = 1
        T[2, 2] = 1
        T[3, 3] = 1
        T[:-1, -1] = ts_f

        # Final affine matrix
        AFF = torch.matmul(T, torch.matmul(Sh, torch.matmul(Sc, torch.matmul(Rz, torch.matmul(Ry, Rx)))))

        # Final vox2vox matrix
        atlas_aff_combined = torch.matmul(AFF, self.aff_A)
        VOX2VOX = torch.matmul(torch.inverse(atlas_aff_combined), self.aff_I)

        # # Nonlinear part: upscale field
        grids_field = self.grids.to(dtype=self.dtype, copy=True)
        grids_field = torch.unsqueeze(grids_field, 0)
        grids_field = grids_field.permute(0, 2, 3, 4, 1)
        for i in range(3):
            grids_field[:, :, :, :, i] = self.FIELD_siz[i]*(grids_field[:, :, :, :, i] / (self.siz_tensor[i]))

        # do not prefilter (field contains the control coefficients, not the displacement values)
        FIELD_rearranged = torch.unsqueeze(field_f.permute(3, 0, 1, 2), dim=0)
        field_resampled = interpol.grid_pull(FIELD_rearranged, grids_field, interpolation=3, prefilter=False, extrapolate=True, bound='dft')
        field_resampled = torch.squeeze(field_resampled)

        # Prepare the new grid for interpolation
        grids_new_atlas = torch.zeros_like(self.grids, dtype=self.dtype)
        for d in range(3):
            grids_new_atlas[d, :, :, :] = VOX2VOX[d, 0] * self.grids[0, :, :, :] \
                                        + VOX2VOX[d, 1] * self.grids[1, :, :, :] \
                                        + VOX2VOX[d, 2] * self.grids[2, :, :, :] \
                                        + VOX2VOX[d, 3] \
                                        + field_resampled[d, :, :, :]

        grids_new_atlas = torch.unsqueeze(grids_new_atlas, 0)
        grids_new_atlas = grids_new_atlas.permute(0, 2, 3, 4, 1)

        return grids_new_atlas, field_f

    def get_grad_penalty(self, field_resampled):

        # For the regularizer: bear in mind that the voxel sizes cancel each other
        # (numerator: because field is in voxels, you need to multiply by voxel size)
        # (denominator: because spacing is also in voxels, you need to multiply by voxel size)
        dy = field_resampled[:, 1:, :, :] - field_resampled[:, :-1, :, :]
        dx = field_resampled[:, :, 1:, :] - field_resampled[:, :, :-1, :]
        dz = field_resampled[:, :, :, 1:] - field_resampled[:, :, :, :-1]

        dy = dy[:, ::self.skip, ::self.skip, ::self.skip]
        dx = dx[:, ::self.skip, ::self.skip, ::self.skip]
        dz = dz[:, ::self.skip, ::self.skip, ::self.skip]

        dy = dy * dy
        dx = dx * dx
        dz = dz * dz
        grad_penalty = (torch.mean(dx) + torch.mean(dy) + torch.mean(dz)) / 3.0

        return grad_penalty

    def get_grad_penalty_new(self, field_voxels):
        """Expects a (nx, ny, nz, 3) tensor, in voxels"""
        # divide by 2 because it's the log of a Gaussian distribution
        return (field_voxels * membrane3(field_voxels)).mean() / 2

    def get_priors(self, atlas_sampling_grid, return_llhood=True):

        atlas_resampled = np.zeros([*self.Iskip.shape, self.n_classes])

        atlas_class = torch.as_tensor(self.A_rearranged[:,0,:,:,:], device=self.device, dtype=self.dtype)
        atlas_class_resampled = interpol.grid_pull(atlas_class, atlas_sampling_grid, interpolation=1)
        atlas_class_resampled = torch.squeeze(atlas_class_resampled).clamp(0, 1)
        if return_llhood:
            LHOODS = torch.zeros([*self.Iskip.shape, self.n_classes], dtype=self.dtype, device=self.device)
            LHOODS[:, :, :, 0] = torch.squeeze(atlas_class_resampled) * self.GAUSSIAN_LHOODS[:, :, :, 0]
        else:
            normalizer = self.GAUSSIAN_LHOODS[:, :, :, 0] * torch.squeeze(atlas_class_resampled) + torch.tensor(1e-6,device=self.device, dtype=self.dtype)

        atlas_class_resampled = atlas_class_resampled.detach()
        atlas_resampled[:,:,:,0] = atlas_class_resampled.cpu().numpy()
        torch.cuda.empty_cache()

        for c in range(self.n_classes-1):
            atlas_class = torch.as_tensor(self.A_rearranged[:,c+1,:,:,:], device=self.device, dtype=self.dtype)
            atlas_class_resampled = interpol.grid_pull(atlas_class, atlas_sampling_grid, interpolation=1)
            atlas_class_resampled = torch.squeeze(atlas_class_resampled).clamp(0, 1)
            num_gaussians = self.gmm_components[c].int()
            for g in range(num_gaussians):
                gaussian_number = sum(self.gmm_components[:c].int()) + g
                if return_llhood:
                    LHOODS[:, :, :, c + 1] += torch.squeeze(atlas_class_resampled) * self.GAUSSIAN_LHOODS[:, :, :, gaussian_number + 1]
                else:
                    normalizer += self.weights[gaussian_number] * self.GAUSSIAN_LHOODS[:, : , :, gaussian_number + 1] * torch.squeeze(atlas_class_resampled)

            atlas_class_resampled = atlas_class_resampled.detach()
            atlas_resampled[:,:,:,c+1] = atlas_class_resampled.cpu().numpy()
            torch.cuda.empty_cache()

        if return_llhood:
            return LHOODS, atlas_resampled
        else:
            return normalizer, atlas_resampled
