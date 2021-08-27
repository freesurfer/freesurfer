function [y,u,s,v] = fast_spatcornoise(Nf,Nv,pvs,v)
% [y,u,s,v] = fast_spatcornoise(Nf,Nv,<pvs>,<v>)
%
% Create spatially correlated, temporally white noise.
%
% Nf = number of frames
% Nv = number of voxels
% pvs = percent variance spanned. If not present or empty,
%   pvs = ones(Nv,1), ie white. pvs will be normalized to 
%   sum to 1. It can be computed from s with pvs = diag(s), 
%   so it can be obtained from a previous call.
% v = spatial eigenvectors. This can be obtained from a 
%   previous call in order to generate a new set of noise
%   with the same spatial correlation. The spatial correlation
%   will be v*s*s*v'
%
% y = (u*s)*v'; 
%
% y is rescaled globally so that the variance over all elements
% is 1.
%
%
%


%
% fast_spatcornoise.m
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

if(nargin < 2 | nargin > 4)
  fprintf('[y,u,s,v] = fast_spatcornoise(Nf,Nv,<pvs>,<v>)\n');
  return;
end

if(~exist('pvs','var')) pvs = []; end
if(isempty(pvs)) 
  % This will make noise spatially white
  pvs = ones(Nv,1);
end

pvs = pvs(:); % column vector
pvs = pvs/sum(pvs);
pvs = sort(pvs);
pvs = flipud(pvs); % largest first

ydim = length(find(pvs>0));
ydim = min(ydim,Nf);
ydim = min(ydim,Nv);

s = diag(pvs(1:ydim));

[u blah blah] = fast_svd(randn(Nf,ydim));
if(~exist('v','var')) 
  [blah blah v] = fast_svd(randn(ydim,Nv));
end

y = (u*s)*v';

% Rescale so that the variance measured over the
% entire data set is 1.
yss2 = std(y(:));
y = y/yss2;

% Rescale s so that y = (u*s)*v' is still true
s = s/yss2;

return;






