function acf = fast_arnw_acf(phi,nf,alpha)
%
% acf = fast_arnw_acf(phi,nf,alpha)
% 
% Computes the autocorrelation function of an AR(N)+White noise
% process given the AR coefficients and alpha.
%
% phi = [phi1 phi2 ... phiN]; 
%   phi1-phiN are the ARN parameters
%   Eg, for an AR1, phi1 = 0.5 corresponds
%   to an acf = (0.5).^[1:nf]
%
% phi has a different form than that resulting from aryule, arburg, 
% arcov, or armcov. These functions return a vector of the form 
% [1 -phi].
%
% If alpha = 0, [], or is not present, then it is
% ignored. Otherwise, the ARN acf is scaled by (1-alpha),
% keeping the 0 lag at 1.
%
% Pure ARN process:
%   y(j) = phi(2)*y(j-1) + phi(3)*y(j-2) ... + e(j)
%
% nf is the number of delays over which to compute
%   the autocorreation function.
%
% Unstable if any abs(pole) is > 1, where poles = roots([1 -phi]);
% 
% See also: fast_arnw_fiterr
%
% Ref: Time Series Analysis by Hamilton, pg. 59.
%
%


%
% fast_arnw_acf.m
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

acf = [];

if(nargin < 2 | nargin > 3)
  fprintf('acf = fast_arnw_acf(phi,nf,alpha)\n');
  return;
end

if(~exist('alpha','var')) alpha = []; end
if(isempty(alpha))        alpha = 0;  end

order = length(phi);
phi = [1 phi(:)'];

F(1,:) = phi(2:order+1); 
F(2:order,:) = [eye(order-1) zeros(order-1,1)];
G = inv(eye(order.^2) - kron(F,F));
g = G(1:order,1);
g2 = g/g(1);

acf = g2;
for j=order+1:nf
  acf(j) = 0;
  for m=1:order
    acf(j) = acf(j) + acf(j-m)*phi(m+1) ;
  end
end

acf = acf(:);

if(alpha > 0)
  acf(2:nf) = (1-alpha)*acf(2:nf);
end

return
