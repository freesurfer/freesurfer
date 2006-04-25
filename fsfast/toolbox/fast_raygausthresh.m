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
% $Id: fast_raygausthresh.m,v 1.1 2006/04/25 01:31:41 greve Exp $

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
