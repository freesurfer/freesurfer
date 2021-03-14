function y = fast_monotonize(x)
% y = fast_monotonize(x)
%
% Makes columns of x monotonically decreasing
%
%


%
% fast_monotonize.m
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

if(nargin ~= 1)
  fprintf('y = fast_monotonize(x)\n');
  return;
end

nf = size(x,1);

y = x;
for n = 2:nf
  ind = find(y(n,:) > y(n-1,:));
  y(n,ind) = y(n-1,ind) ;
end

return;
