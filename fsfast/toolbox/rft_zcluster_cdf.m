function cdf = rft_zcluster_cdf(csize,vthresh,fwhm,ssize,D)
% cdf = rft_zcluster_cdf(csize,vthresh,fwhm,ssize,D)
%
% Prob that a cluster >= csize defined by threshold vthresh will be
% found in a D-dim z-field with fwhm smoothness of ssize search space.
% vthresh is the voxel-wise z thresh


%
% rft_zcluster_cdf.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
%    $Revision: 1.3 $
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

cdf = [];
if(nargin ~= 5)
  fprintf('cdf = rft_zcluster_cdf(csize,vthresh,fwhm,ssize,D)\n');
  return;
end

W = fwhm/sqrt(4*log(2));
dLh = W.^(-D); % Same as FSL?
%dLh = 0.493971;

if(0)
Em = ssize .* ((2*pi).^(-(D+1)/2)) .* dLh .* (vthresh.^(D-1)) .* ...
     exp(-(vthresh.^2)/2); 
else
% This appears to be the FSL equation from the techrep
Em = ssize .* ((2*pi).^(-(D+1)/2)) .* dLh .* (vthresh.^(D-1) - 1) .* ...
     exp(-(vthresh.^2)/2); 
end

beta = ( ((gamma(D/2+1) .* Em)) ./ (ssize.*mri_zcdf(-vthresh))).^(2./D);

% Prob than n >= k, ie, the number of voxels in a cluster >= csize
Pnk = exp(-beta.*(csize.^(2./D)));
cdf = 1 - exp(-Em.*Pnk);

%keyboard

return;


