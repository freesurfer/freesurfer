function [p,pR,pG] = fast_raygaus(y,alpha,Rmu,Gmu,Gstd)
% [p pR pG] = raygaus(y,alpha,Rmu,Gmu,Gstd)
% [p pR pG] = raygaus(y,params)
%
% Computes the probability of getting y from
% a Rayleigh-Gaussian mixture.
%
% alpha - fraction in rayleigh
% Rmu - mean of rayleigh
% Gmu - mean of gaussian
% Gstd - std of gaussian
%
% params = [alpha Rmu Gmu Gstd];
%
% p = alpha*pR + (1-alpha)*pG
%
% See also: pdf_rayleigh, gaussian, fast_raygausinit, 
%   fast_raygauscost.
%
%


%
% fast_raygaus.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:05 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

p = [];
if(nargin ~= 2 & nargin ~= 5)
  fprintf('[p pR pG] = fast_raygaus(y,alpha,Rmu,Gmu,Gstd)\n');
  fprintf('[p pR pG] = fast_raygaus(y,params)\n');
  return;
end

if(nargin == 2)
  params = alpha;
  alpha = params(1);
  Rmu   = params(2);
  Gmu   = params(3);
  Gstd   = params(4);
end

pR = pdf_rayleigh(y,Rmu);
pG = gaussian(y,Gmu,Gstd);

p = alpha*pR + (1-alpha)*pG;

return;
