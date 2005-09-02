function [volsm, G] = fast_smooth3d(vol,cfwhm,rfwhm,sfwhm)
% [volsm, G] = fast_smooth3d(vol,cfwhm,rfwhm,sfwhm)
% 
% 3D gaussian smoother.
%
% cfwhm - fwhm for cols
% rfwhm - fwhm for rows
% sfwhm - fwhm for slice
%
% Note: does not attempt to handle wrap-around
%
% $Id: fast_smooth3d.m,v 1.4 2005/09/02 17:30:47 greve Exp $

if(nargin ~= 4)
  fprintf('[volsm, G] = fast_smooth3d(vol,cfwhm,rfwhm,sfwhm)\n');
  return;
end

if(isfield(vol,'vol')) vol = vol.vol; end

rstd = rfwhm/sqrt(log(256.0));
cstd = cfwhm/sqrt(log(256.0));
sstd = sfwhm/sqrt(log(256.0));

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

volsm = vol;

return
