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
% $Id: arCorrFun.m,v 1.1 2003/03/04 20:47:33 greve Exp $

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

