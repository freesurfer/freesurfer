function [r, mask, voxmain, voxghost] = fast_ghostdet(vol, mask)
% [r mask voxmain voxghost] = fast_ghostdet(vol, <mask>)
%
% vol is [nrows ncols nslices nframes]
%   Each slice is an epi slice
%   The phase encode changes from one row to another (beware
%     matlabs row major!)
% mask is a binary mask. If unspecified, a mask is created based
%   on the number of voxels above threshold.
%
% r is the ratio of the mean of the voxels in the main to the
% mean of the voxels in the ghost (excluding voxels that are in both).
%
% voxmain and voxghost are vectors of voxels in the main image
% and their corresponding voxels in the ghost (again excluding
% voxels that are in both).
% 
%


%
% fast_ghostdet.m
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

% To test: this should force r=1
if(0)
  mask = zeros(64,64,5);
  mask(15:50,15:50,:) = 1;
  vol = randn(64,64,5,20);
  volfft = fft2(vol);
  volfft(2:2:end,:,:,:) = 0;
  vol = abs(ifft2(volfft));
end

r = [];

if(nargin < 1 | nargin > 2)
  fprintf('[r, mask] = fast_ghostdet(vol, <mask>)\n');
  return;
end

[nrows ncols nslices nframes] = size(vol);
nrows2 = round(nrows/2); % half the number of rows
nvoxels = nrows*ncols*nslices;

if(exist('mask','var'))
  % Mask passed as an argument, check size
  if(size(mask,1) ~= nrows | size(mask,2) ~= ncols | ...
     size(mask,3) ~= nslices)
    fprintf('ERROR: fast_ghostdet: mask size mismatch\n');
    return;
  end
else
  % Create Mask
  volmean = mean(vol,4); % mean over frames
  globalmean = mean(reshape(volmean,[nvoxels 1])); % mean over everything
  thresh = 0.5*globalmean;
  mask = volmean > thresh;
  % Should probably dilate mask
end

% Get a mask of the ghost
ghostmask = fast_ghostmask(mask);

% Get a mask of the voxels in the mask but not in the ghost
masknotghost = mask & ~ghostmask;
if(length(find(masknotghost)) == 0)
  fprintf('ERROR: fast_ghostdet: cannot separate main from mask\n');
  return;
end

% Now Get a mask of the voxels in the ghost but not in the mask
% ghostnotmask = fast_ghostmask(masknotghost);
% This is actually not helpful.

% Reshape vol to a more convenient form
vol = reshape(vol,[nrows*ncols nslices nframes]);

% Go through each slice, accumulating voxels
voxmain  = [];
voxghost = [];
for slice = 1:nslices

  % linear indices in slice
  ind1 = find(masknotghost(:,:,slice));
  nhits = length(ind1);

  % r1 and c1 are the rows and corresponding columns of the
  % voxels in the mask but not in the ghost.
  [r1 c1] = ind2sub(size(mask),ind1);

  % Now need to compute where these voxels go to in the ghost. Only
  % the row index changes, the columns stay the same.
  indr1TH = find(r1 <= nrows2); % top half
  indr1BH = find(r1 > nrows2);  % bottom half

  % Compute the ghost row numbers
  r2 = zeros(size(r1));
  r2(indr1TH) = r1(indr1TH) + nrows2;
  r2(indr1BH) = r1(indr1BH) - nrows2;

  % Convert the ghost row/col indices to linear indices into the slice
  ind2 = sub2ind(size(mask),r2,c1);

  voxmainslice = squeeze(vol(ind1,slice,:));
  voxmainslice = reshape(voxmainslice,[nhits*nframes 1]);
  voxmain = [voxmain; voxmainslice];
  
  voxghostslice = squeeze(vol(ind2,slice,:));
  voxghostslice = reshape(voxghostslice,[nhits*nframes 1]);
  voxghost = [voxghost; voxghostslice];
  
end

mainmean  = mean(voxmain);
ghostmean = mean(voxghost);

r = mainmean/ghostmean;

return;



