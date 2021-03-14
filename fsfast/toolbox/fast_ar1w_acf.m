function acf = fast_ar1w_acf(alpha, rho, nf)
% acf = fast_ar1w_acf(alpha, rho, nf)
%
% Autocorrelation function of AR1+White Noise
%   acf(0) = 1
%   acf(n) = (1-alpha)*rho^n


%
% fast_ar1w_acf.m
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

if(nargin ~= 3)
  fprintf('acf = fast_ar1w_acf(alpha, rho, nf)\n');
  return;quit;
end

nv = size(alpha,2);

nn = repmat([0:nf-1]',[1 nv]);
rho = repmat(rho,[nf 1]);
acf = (rho.^nn) .* repmat((1-alpha),[nf 1]);
acf(1,:) = 1;

return;
