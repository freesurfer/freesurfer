function m = reshape2d(x)
% m = reshape2d(x)
% Reshape into a 2d array where
% size(m,1) = size(x,1)


%
% reshape2d.m
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
  msg = 'm = reshape2d(x)'
  qoe(msg);error(msg);
end

szx = size(x);
xdim = length(szx);

if(xdim < 2)
  msg = sprintf('Dimension of x is only %d',xdim);
  qoe(msg);error(msg);
end

m = reshape(x, [szx(1) prod(szx(2:xdim))]);
return
