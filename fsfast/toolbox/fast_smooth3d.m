function volsm = fast_smooth3d(vol,cfwhm,rfwhm,sfwhm,UseBB)
% volsm = fast_smooth3d(vol,cfwhm,rfwhm,sfwhm,UseBB)
% 
% 3D gaussian smoother.
%
% cfwhm - fwhm for cols
% rfwhm - fwhm for rows
% sfwhm - fwhm for slice
% UseBB - reduce size to bounding box of non-zero voxels padded
% with 4STDs in each dim. This will not yield exactly the same 
% result, but it will be close, but it can speed things up a lot.

%
% fast_smooth3d.m
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

volsm = [];

if(nargin ~= 4 & nargin ~= 5)
  fprintf('volsm = fast_smooth3d(vol,cfwhm,rfwhm,sfwhm,<UseBB>)\n');
  return;
end

if(~exist('UseBB','var')) UseBB = []; end
if(isempty(UseBB))        UseBB = 0; end

if(isfield(vol,'vol')) vol = vol.vol; end

rstd = rfwhm/sqrt(log(256.0));
cstd = cfwhm/sqrt(log(256.0));
sstd = sfwhm/sqrt(log(256.0));

if(UseBB)
  % This is not exactly the same near the edges of the volume
  [nrows0 ncols0 nslices0 nframes0] = size(vol);
  nvox0 = prod(size(vol));
  bb = mri_boundingbox(vol);
  %fprintf('bb min %d %d %d max %d %d %d\n',bb(1,:),bb(2,:));
  nstds = 4;
  rmin = max(round(bb(1,1)-nstds*rstd),1);
  rmax = min(round(bb(2,1)+nstds*rstd),nrows0);
  cmin = max(round(bb(1,2)-nstds*cstd),1);
  cmax = min(round(bb(2,2)+nstds*cstd),ncols0);
  smin = max(round(bb(1,3)-nstds*sstd),1);
  smax = min(round(bb(2,3)+nstds*sstd),nslices0);
  vol = vol(rmin:rmax,cmin:cmax,smin:smax,:);
  nvox = prod(size(vol));
  %fprintf('nvox0/nvox = %g\n',nvox0/nvox);
end

[nrows ncols nslices nframes] = size(vol);

% Do the rows:
if(rstd > 0)
  vol = reshape(vol,[nrows ncols*nslices*nframes]);
  vol = fast_smooth1d(vol,rstd);
  vol = reshape(vol,[nrows ncols nslices nframes]);
end

% Do the cols:
if(cstd > 0)
  vol = permute(vol,[2 1 3 4]);
  vol = reshape(vol,[ncols nrows*nslices*nframes]);
  vol = fast_smooth1d(vol,cstd);
  vol = reshape(vol,[ncols nrows nslices nframes]);
  vol = permute(vol,[2 1 3 4]);
end

% Do the slices
if(sstd > 0)
  vol = permute(vol,[3 2 1 4]);
  vol = reshape(vol,[nslices ncols*nrows*nframes]);
  vol = fast_smooth1d(vol,sstd);
  vol = reshape(vol,[nslices ncols nrows nframes]);
  vol = permute(vol,[3 2 1 4]);
end

if(~UseBB)
  volsm = vol;
else
  volsm = zeros(nrows0,ncols0,nslices0,nframes0);
  volsm(rmin:rmax,cmin:cmax,smin:smax,:) = vol;
end


return
