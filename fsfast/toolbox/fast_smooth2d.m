function [volsm, G] = fast_smooth2d(vol,cfwhm,rfwhm,sfwhm)
% [volsm, G] = fast_smooth2d(vol,cfwhm,rfwhm,sfwhm)
% 
% 3D gaussian smoother.
%
% cfwhm - fwhm for cols
% rfwhm - fwhm for rows
% sfwhm - fwhm for slice
%
% Note: does not attempt to handle wrap-around
%
%


%
% fast_smooth2d.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

if(nargin ~= 3)
  fprintf('[volsm, G] = fast_smooth2d(vol,cfwhm,rfwhm)\n');
  return;
end

if(isfield(vol,'vol')) vol = vol.vol; end

rstd = rfwhm/sqrt(log(256.0));
cstd = cfwhm/sqrt(log(256.0));

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

volsm = vol;

return
