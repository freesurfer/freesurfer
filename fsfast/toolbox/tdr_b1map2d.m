function [b1map, mask] = tdr_b1map2d(body,head)
% [b1map mask] = tdr_b1map2d(body,head)
% body = complex nrows x ncols image from body coil
% head = complex nrows x ncols image from single head coil
%
% Issues:
%  1. Should probably make mask from head, not body
%  2. Should scale be changed in some way?
%
%


%
% tdr_b1map2d.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

rthresh = 1; % relative threshold for masking
nsmooth = 2; % number of in-plane voxels to smooth over
nfit = 20; % number of in-mask voxels to use to fit edge voxels

b1map = [];
mask = [];
if(nargin ~= 2)
  fprintf('[b1map mask] = tdr_b1map2d(body,head)\n');
  return;
end

img = abs(body);
szimg   = size(img);
nrows   = szimg(1);
ncols   = szimg(2);
nvox = prod(szimg);

imgmn = mean(img(:));
mask = img > rthresh*imgmn;
indmask = find(mask);
nmask = length(indmask);
fprintf('nmask = %d  %g\n',nmask,nmask/nvox);

% Smooth mask in 2d
mask_smth = fast_smooth2d(mask,nsmooth,nsmooth);

% Compute complex Head to Body ratio
rhb = head./body;

% We want to smooth the rhb, but only in the mask

% Masked Head to Body Ratio
rhbm = rhb .* mask;

% Smoothed Masked Head to Body Ratio
rhbm_smth = fast_smooth2d(rhbm,nsmooth,nsmooth);

% Rescale to account for mask edge effects in smoothing
rhbm_smth_sc = rhbm_smth;
rhbm_smth_sc(indmask) = rhbm_smth(indmask)./mask_smth(indmask);

b1map0 = rhbm_smth_sc;

% Now need to extrapolate to voxels near the edge of mask
[rmask cmask] = ind2sub(size(mask),indmask);

% Dilate mask by 3 voxels in-plane. This defines the edge voxels
m3d = fast_dilate(mask,3,0,1);
medge = m3d & ~mask;
indedge = find(medge);
nedge = length(indedge);
fprintf('nedge = %d  %g\n',nedge,nedge/nvox);
[redge cedge] = ind2sub(size(medge),indedge);

% For each voxel in the edge, find the nfit closest voxels that are
% also in the mask. Then do a linear fit to the b1map in the mask,
% and then extrapolate the fit the edge voxel.
b1map = b1map0;
for n=1:nedge
  d2 = (redge(n)-rmask).^2 + (cedge(n)-cmask).^2;
  [d2sort isort] = sort(d2);
  rfit = rmask(isort(1:nfit));
  cfit = cmask(isort(1:nfit));
  X = [ones(nfit,1) rfit cfit];
  y = b1map0(indmask(isort(1:nfit)));
  beta = inv(X'*X)*X'*y;
  yn = [1 redge(n) cedge(n)]*beta;
  b1map(indedge(n)) = yn;
end

return;

