function eulerangles = MRImdc2euler(Mdc)
% eulerangles = MRImdc2euler(Mdc)
%
% Solves for the Euler angles given the 3x3 matrix of direction cosines.
% This code does not work in all cases.
%


%
% MRImdc2euler.m
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


eulerangles = [];
if(nargin ~= 1)
  fprintf('eulerangles = MRImdc2euler(Mdc)\n');
  return;
end

theta = acos(Mdc(3,3));
phi = asin(Mdc(3,2)/(sin(theta)+eps));
psi = asin(Mdc(2,3)/(sin(theta)+eps));

eulerangles = [theta phi psi];

return;



