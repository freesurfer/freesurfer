#=
  Original Author: Anastasia Yendiki

  Copyright © 2022 The General Hospital Corporation (Boston, MA) "MGH"
 
  Terms and conditions for use, reproduction, distribution and contribution
  are found in the 'FreeSurfer Software License Agreement' contained
  in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 
  https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 
  Reporting: freesurfer@nmr.mgh.harvard.edu
=#

using LinearAlgebra, Statistics

export DTI, dti_fit, dti_write


"Container for outputs of a DTI fit"
struct DTI
  eigval1::MRI
  eigval2::MRI
  eigval3::MRI
  eigvec1::MRI
  eigvec2::MRI
  eigvec3::MRI
  rd::MRI
  md::MRI
  fa::MRI
end


"""
    dti_fit_ls(dwi::MRI, mask:MRI)

Fit tensors to DWIs and return a `DTI` structure.
"""
function dti_fit(dwi::MRI, mask::MRI)
  dti_fit_ls(dwi::MRI, mask::MRI)
end


"""
    dti_fit_ls(dwi::MRI, mask:MRI)

Perform least-squares fitting of tensors from DWIs and return a `DTI` structure.
"""
function dti_fit_ls(dwi::MRI, mask::MRI)

  if isempty(dwi.bval)
    error("Missing b-value table from input DWI structure")
  end

  if isempty(dwi.bvec)
    error("Missing gradient table from input DWI structure")
  end

  ib0 = (dwi.bval .== minimum(dwi.bval))

  A = hcat(dwi.bvec[:,1].^2, 2*dwi.bvec[:,1].*dwi.bvec[:,2],
           2*dwi.bvec[:,1].*dwi.bvec[:,3], dwi.bvec[:,2].^2,
           2*dwi.bvec[:,2].*dwi.bvec[:,3], dwi.bvec[:,3].^2)
  A = hcat(-dwi.bval .* A, ones(size(A, 1), 1))

  pA = pinv(A)

  Eval1 = MRI(mask, 1)
  Eval2 = MRI(mask, 1)
  Eval3 = MRI(mask, 1)
  Evec1 = MRI(mask, 3)
  Evec2 = MRI(mask, 3)
  Evec3 = MRI(mask, 3)

  Threads.@threads for iz in 1:size(dwi.vol, 3)
    for iy in 1:size(dwi.vol, 2)
      for ix in 1:size(dwi.vol, 1)
        mask.vol[ix, iy, iz] == 0 && continue

        # Only use positive DWI values to fit the model
        ipos = dwi.vol[ix, iy, iz, :] .> 0
        npos = sum(ipos)

        if npos == length(dwi.bval)
          D = pA * log.(dwi.vol[ix, iy, iz, :])
        elseif npos > 6
          sum(ipos .&& ib0) == 0 && continue
          D = pinv(A[ipos, :]) * log.(dwi.vol[ix, iy, iz, ipos])
        else
          continue
        end

        E = eigen([D[1] D[2] D[3];
                   D[2] D[4] D[5];
                   D[3] D[5] D[6]])

        Eval1.vol[ix, iy, iz]    = E.values[3]
        Eval2.vol[ix, iy, iz]    = E.values[2]
        Eval3.vol[ix, iy, iz]    = E.values[1]
        Evec1.vol[ix, iy, iz, :] = E.vectors[:, 3]
        Evec2.vol[ix, iy, iz, :] = E.vectors[:, 2]
        Evec3.vol[ix, iy, iz, :] = E.vectors[:, 1]
      end
    end
  end

  return DTI(Eval1, Eval2, Eval3, Evec1, Evec2, Evec3,
             dti_maps(Eval1, Eval2, Eval3)...)
end

"""
    dti_maps(eigval1::MRI, eigval2::MRI, eigval3::MRI)

Compute radial diffusivity (RD), mean diffusivity (MD), and fractional
anisotropy (FA) maps from the 3 eigenvalues the diffusion tensors.

Return RD, MD, and FA maps as MRI structures.
"""
function dti_maps(eigval1::MRI, eigval2::MRI, eigval3::MRI)

  rd = MRI(eigval1)
  md = MRI(eigval1)
  fa = MRI(eigval1)

  imask = (eigval1.vol .!= 0)

  rd.vol[imask] = eigval1.vol[imask] + eigval2.vol[imask]
  md.vol[imask] = (rd.vol[imask] + eigval3.vol[imask]) / 3
  rd.vol[imask] /= 2

  fa.vol[imask] = sqrt.(((eigval1.vol[imask] - md.vol[imask]).^2 +
                         (eigval2.vol[imask] - md.vol[imask]).^2 +
                         (eigval3.vol[imask] - md.vol[imask]).^2) ./
                         (eigval1.vol[imask].^2 +
                          eigval2.vol[imask].^2 +
                          eigval3.vol[imask].^2) * 3/2)

  return rd, md, fa
end


"""
    dti_write(dti::DTI, basename::String)

Write the volumes from a `DTI` structure that was created by `dti_fit()`
to files whose names start with the specified base name.
"""
function dti_write(dti::DTI, basename::String)

  for var in fieldnames(DTI)
    mri_write(getfield(dti, var), basename * "_" * string(var) * ".nii.gz")
  end
end


