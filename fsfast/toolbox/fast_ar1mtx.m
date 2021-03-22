function M = fast_ar1mtx(rho,len)
% M = fast_ar1mtx(rho,len)
%
% Produces a toeplitz AR1 matrix
%


%
% fast_ar1mtx.m
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

M = [];

if(nargin ~= 2)
  fprintf('USAGE: M = fast_ar1mtx(rho,len)\n');
  return;
end

M = [];
for c = 0:len-1;
  ind = abs([-c : len-c-1]);
  r = rho .^ ind;
  M = [M; r];
end

return;
