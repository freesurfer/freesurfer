function R = arCorrFun(alpha, rho, nMaxDelay, Sided)
%
% function R = arCorrFun(alpha, rho, nMaxDelay, <Sided>)
%
% nMaxDelay + 1 = number of components in the One-Sided ACor.
%
% The alpha-rho model of the correlation function:
% r(0) = 1
% r(n) = (1-alpha)*rho^n,  0 < n <= nMaxDelay
%
% Defaults: Sided = 1
%
%


%
% arCorrFun.m
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

if(nargin == 3)   Sided = 1; end % default

nR = nMaxDelay + 1; 

if(Sided == 1)
  R=zeros(nR,1);
  R(1)=1;
  R(2:nR) = (1-alpha)*(rho.^([1:nR-1]));
else
  R=zeros(2*(nR-1)+1,1);
  R(nR)=1;
  x = [1:nR-1];
  z = (1-alpha)*(rho.^(x));
  R(x) = z(nR-x);
  R(nR+x) = z(x);
end

return;

