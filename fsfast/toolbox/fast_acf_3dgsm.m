function [acf, kret] = fast_acf_3dgsm(fwhm,lag,ndim)
% 
% [acf, kret] = fast_acf_3dgsm(fwhm,<lag>,<ndim>)
%
% Computes the expected spatial autocorrelation function at lag lag
% for white noise smoothed isotropically with a 3D gaussian filter
% fwhm voxels. Note: gauss std = fwhm/2.36.
%
% If lag is not specified or null, then lag = 1.
% lag can be a list of integers.
%
% ndim is the dimension of the gaussian smoothing. Default ndim=3
% I think that for isotropic smoothing, ndim=1 should give the same
% answer as ndim=3. For this reason, I did not implement ndim=2.
%
%
%
% (c) Douglas N. Greve, 2004.


%
% fast_acf_3dgsm.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
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

acf = [];
kret = [];

if(nargin < 1 | nargin > 3)
  fprintf('ar1 = ar1gaussian(fwhm,<lag>)\n');
  return;
end

if(~exist('lag','var')) lag = []; end
if(isempty(lag)) lag = 1; end

if(~exist('ndim','var')) ndim = []; end
if(isempty(ndim)) ndim = 3; end

gstd = fwhm/2.36;
w = round(4*fwhm) + max(lag);
x = -w:w;
nx = length(x);

k1 = gaussian(x,0,gstd);
if(ndim == 1)
  k1 = k1/sum(k1);
  ar0 =  sum(reshape1d((k1.*k1)));
  nth = 1;
  for n = lag
    acf(nth) = sum(reshape1d((k1(1:end-n)).*(k1(n+1:end))));
    nth = nth + 1;
  end
  acf = acf/ar0;
  kret = k1;
  return;
end


k2 = k1' * k1;

k3 = reshape(reshape1d(k2) * k1,[nx nx nx]);
k3 = k3/sum(reshape1d(k3));

ar0 =  sum(reshape1d((k3.*k3)));

nth = 1;
for n = lag
  acf(nth) = sum(reshape1d((k3(1:end-n,:,:)).*(k3(n+1:end,:,:))));
  nth = nth + 1;
end

acf = acf/ar0;

return;








