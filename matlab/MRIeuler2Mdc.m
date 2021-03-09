function Mdc = MRIeuler2Mdc(eulerangles)
% Mdc = MRIeuler2Mdc(eulerangles)
% 
% Returns the matrix of direction cosines given the Euler angles.
% Not sure if this works.
%


%
% MRIeuler2Mdc.m
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


Mdc = [];
if(nargin ~= 1)
  fprintf('Mdc = MRIeuler2Mdc(eulerangles)\n');
  return;
end

theta = eulerangles(1);
phi   = eulerangles(2);
psi   = eulerangles(3);

Mdc(1,1) = cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi);
Mdc(2,1) = sin(phi)*cos(psi) + cos(phi)*cos(theta)*sin(psi);
Mdc(3,1) = -cos(phi)*sin(theta);

Mdc(1,2) = -sin(phi)*cos(theta)*cos(psi) - cos(phi)*sin(psi);
Mdc(2,2) = -sin(phi)*cos(theta)*sin(psi) + cos(phi)*cos(psi);
Mdc(3,2) = sin(theta)*sin(phi);

Mdc(3,1) = sin(theta)*cos(psi);
Mdc(3,2) = sin(theta)*sin(psi);
Mdc(3,3) = cos(theta);

return;
