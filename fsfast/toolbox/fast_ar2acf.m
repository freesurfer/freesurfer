function acf = fast_ar2acf(phi,nmax)
%
% Computes the autocorrelation function given the AR 
% coefficients (phi)
%
% acf = fast_acf_yw(phi)
%
% phi is the result of aryule, arburg, arcov, or armcov. 
%   phi length is order+1, and phi(1) is always 1.
%
% y(j) = -phi(2)*y(j-1) - phi(3)*y(j-2) ... + e(j)
%
% nmax is the number of delays over which to compute
%   the autocorreation function.
%
% Unstable if abs of abs(pole) is > 1
% poles = roots(phi);
% 
% Ref: Time Series Analysis by Hamilton, pg. 59.
%
%


%
% fast_ar2acf.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
%    $Revision: 1.5 $
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

acf = [];

if(nargin ~= 2)
  fprintf('USAGE: acf = fast_acf_yw(phi,nmax)\n');
  return;
end

order = length(phi) - 1;

% Negative is needed here %
phi = -phi;
phi(1) = 1;

F(1,:) = phi(2:order+1); 
F(2:order,:) = [eye(order-1) zeros(order-1,1)];
G = inv(eye(order.^2) - kron(F,F));
g = G(1:order,1);
g2 = g/g(1);

acf = g2;
for j=order+1:nmax
  acf(j) = 0;
  for m=1:order
    acf(j) = acf(j) + acf(j-m)*phi(m+1) ;
  end
end

acf = reshape1d(acf);

return
