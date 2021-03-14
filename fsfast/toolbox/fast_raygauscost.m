function cost = fast_raygauscost(params,y)
% cost = fast_raygauscost(p,y)
% params = [alpha Rmu Gmu Gstd];
%
% Cost function for Rayleigh-Gaussian Mixture Model
%
% [p pR pG] = fast_raygaus(y,alpha,Rmu,Gmu,Gstd);
% cost = -sum(log10(p));
%
% optparams = fminsearch('fast_raygauscost',params,[],y);
%
% See also fast_raygaus, fast_raygauscost.
%
%


%
% fast_raygauscost.m
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

alpha = params(1);
Rmu   = params(2);
Gmu   = params(3);
Gstd   = params(4);

if(alpha < 0 | alpha > 1 | Rmu < 0 | Gmu < 0 | Gstd < 0)
  cost = 10e10;
  return;
end

p = fast_raygaus(y,alpha,Rmu,Gmu,Gstd);
ind = find(p == 0);
if(~isempty(ind))
  cost = 10e10;
  return;
end

cost = -sum(log10(p));

return;










