function ythresh = fast_raygausthresh(alpha,Rmu,Gmu,Gstd)
% ythresh = fast_raygausthresh(alpha,Rmu,Gmu,Gstd)
% ythresh = fast_raygausthresh(params)
%
% Computes the threshold to separate the Rayleigh voxels
% from the Gaussian voxels in a Rayleigh-Gaussian Mixture
% model.
%
% alpha - fraction in rayleigh
% Rmu - mean of rayleigh
% Gmu - mean of gaussian
% Gstd - std of gaussian
%
% params = [alpha Rmu Gmu Gstd];
%
% See also: pdf_rayleigh, gaussian, fast_raygaus, 
%   fast_raygausinit, fast_raygauscost.
%
%


%
% fast_raygausthresh.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
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

p = [];
if(nargin ~= 1 & nargin ~= 4)
  fprintf('ythresh = fast_raygausthresh(alpha,Rmu,Gmu,Gstd)\n');
  fprintf('ythresh = fast_raygausthresh(params)\n');
  return;
end

if(nargin == 1)
  params = alpha;
  alpha = params(1);
  Rmu   = params(2);
  Gmu   = params(3);
  Gstd  = params(4);
end

ymin = .5*Rmu;
ymax = Gmu+3*Gstd;
dy = (ymax-ymin)/1000;
y = [ymin:dy:ymax];

[p pR pG] = fast_raygaus(y,alpha,Rmu,Gmu,Gstd);
pR = alpha*pR;
pG = (1-alpha)*pG;
d = pG-pR;
[mm mi] = min(abs(d));
ythresh = y(mi);

return;
