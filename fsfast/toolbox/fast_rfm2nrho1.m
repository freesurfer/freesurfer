function [beta, rrho1, nrho1, nrho1hat] = fast_rfm2nrho1(R,nrho1)
% [beta rrho1 nrho1 nrho1hat] = fast_rfm2nrho1(R,<nrho1>)
%
% Computes the correction factor to an AR1 value for the given
% residual forming matrix R. Works by synthesizing a noise
% covariance matrix for a series of true AR1 values (nrho1).
% If nrho1 is not supplied, then uses nrho1 = [-.5:.25:.5].
%
% rrho1 are the expected AR1 values when computed from the 
% residual for each of the true AR1 values nrho1.
%
% beta are the two best-fit linear components relating rrho1 and
% nrho1, ie, nrho1hat = [1 rrho1]*beta, or equivalently nrho1hat =
% beta(1) + rho1*beta(2). This is the way to correct for the bias
% induced by estimating the AR1 from the resiual.
%
% This is sort of an empirical version of Keith's method (but a
% little more accurate). The relationship between rrho1 and nrho1
% is pretty linear over that range, at least for event-related
% designs.
% 
%


%
% fast_rfm2nrho1.m
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

beta = [];
if(nargin ~= 1)
  fprintf('[beta rrho1 nrho1 nrho1hat] = fast_rfm2nrho1(R,<nrho1>)\n');
  return;
end

ntp = size(R,1);
nn = [0:ntp-1];

if(~exist('nrho1','var'))  nrho1 = [-.5:.25:.5]'; end
nrho1 = nrho1(:);

nthrho = 1;
tic;
for rho = nrho1'
  acf = rho.^nn;
  M = toeplitz(acf);
  D = R*M*R;
  rrho1(nthrho,1) = mean(diag(D,1));
  %nrho1kjw = fast_yacf_kjw([1 rrho1(nthrho,1)],R);
  fprintf('%2d  %g  %g    (t=%g)\n',nthrho,rho,rrho1(nthrho),toc);
  nthrho = nthrho + 1;  
end

nlist = length(nrho1);

X = [ones(nlist,1) rrho1];
beta = inv(X'*X)*X'*nrho1;
nrho1hat = X*beta;

return;





