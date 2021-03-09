function [wvar, ivar] = fast_flip2wvar(mn1,mn2,var1,var2)
% [wvar, ivar] = fast_flip2wvar(mn1,mn2,var1,var2)
%
% Compute white and instability noise given two measures
% with different proportions. The instability variance is that
% from the first scan.
%

%
% fast_flip2var
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

wvar = [];
ivar = [];

if(nargin ~= 4)
  fprintf('[wvar, ivar] = fast_flip2wvar(mn1,mn2,var1,var2)\n');
  return;
end

% This is the factor by which the instability should change due to
% changing the acq parameters:
v = (mn1./mn2).^2;

wvar = (v.*var2 - var1)./(v-1);
ivar = var1 - wvar;

return;





