function [roivar, mask, navg] = fast_weiskopf_perm(fslice,nroi,mask,navg)
% [roivar, mask, navg] = fast_weiskopf(fslice,nroi,<mask>,<navg>)
%
% Computes the variance over an ROI of a given number of voxels.
%
% fslice - functional data [nr nc nf] = size(fslice);
% nroi - list of number of voxels in each roi to test
% mask - only choose voxels from the given mask. If mask is not
%   specified, then fslice is thresholded to give a mask. If mask
%   is null, then no mask is used.
% navg - randomly select navg sets for each roi, and average the
%   variances across them.
%
% When fslice is white, then the following two curves should be the
% same: loglog(nroi,roivar,nroi,1./nroi)
%
% Results: never indicates that there's a problem because the
% voxels are too dispersed across the volume.
%
%
%
% (c) Douglas N. Greve, 2004


%
% fast_weiskopf_perm.m
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

roivar = [];
polyorder = 2; % detrending order

if(nargin < 2 | nargin > 4)
  fprintf('[roivar, mask, navg] = fast_weiskopf_perm(fslice,nroi,<mask>,<navg>)\n');
  return;
end

[nr nc nf] = size(fslice);
if(nf == 1)
  fprintf('ERROR: fslice only have one frame\n');
  return;
end

if(~exist('navg','var')) navg = []; end
if(isempty(navg)) navg = 10; end

if(~exist('mask','var'))
  % Mask not specified; get by thresholding
  fslicemn = mean(fslice,3);
  fslicegmn = mean(reshape1d(fslicemn));
  mask = fslicemn > 2*fslicegmn;
end
if(isempty(mask)) 
  % Mask specified but null; no mask.
  mask = ones(nr,nc);
end
indmask = find(mask);
nmask = length(indmask);

if(max(nroi) > nmask)
  fprintf('ERROR: max(nroi) > nmask)\n');
  return;
end

% Reshape and mask
fslice = reshape(fslice,[nr*nc nf])';
fslice = fslice(:,indmask);

% Detrend
X = fast_polytrendmtx(1,nf,1,polyorder);
R = eye(nf) - X*inv(X'*X)*X';
fslice = R*fslice;

nrois = length(nroi);
for nthroi = 1:nrois
  nroiuse = nroi(nthroi);
  tmpstd = zeros(navg,1);
  for nthavg = 1:navg
    ind = randperm(nmask);
    ind = ind(1:nroiuse);
    tmpstd(nthavg) = std(mean(fslice(:,ind),2));
  end

  roivar(nthroi) = mean(tmpstd.^2);
  
end

return;

