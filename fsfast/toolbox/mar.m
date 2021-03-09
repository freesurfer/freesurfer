function [m, i] = mar(x)
%
% [m i] = mar(x)
%
% [m i] = max(abs(reshape1d(x)))
%


%
% mar.m
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
  msg = 'USAGE: [m i] = mar(x)';
  qoe(msg);error(msg);
end

n = prod(size(x));

[m i] = max(abs(reshape(x,[n 1])));

return;
