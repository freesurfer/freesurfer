function acf = fast_arnw_acf(p,nf)
%
% acf = fast_arnw_acf(p,nf)
% 
% Computes the autocorrelation function of an AR(N)+White noise
% process given the AR coefficients and alpha.
%
% p = [alpha phi1 phi2 ... phiN]; 
%   phi1-phiN are the ARN parameters consistent.
%   Note: these are negative relative to the "standard" 
%   interpretation. Eg, for an AR1, phi1 = -0.5 corresponds
%   to an acf = (0.5).^[1:nf]
%
% phi is of the form resulting from aryule, arburg, arcov, or armcov. 
%
% Pure ARN process:
%   y(j) = -phi(2)*y(j-1) - phi(3)*y(j-2) ... + e(j)
%
% nf is the number of delays over which to compute
%   the autocorreation function.
%
% Unstable if abs of abs(pole) is > 1
%   poles = roots(phi);
% 
% See also: fast_ar2acf, fast_arnw_fiterr
%
% $Id: fast_arnw_acf.m,v 1.1 2004/05/24 00:08:55 greve Exp $

acf = [];

if(nargin ~= 2)
  fprintf('acf = fast_arnw_acf(p,nf)\n');
  return;
end

alpha = p(1);
phi = p(2:end);
phi = [1 phi(:)'];

acf = fast_ar2acf(phi,nf);
acf(2:end) = (1-alpha)*acf(2:end);


return
