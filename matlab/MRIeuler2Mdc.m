function Mdc = MRIeuler2Mdc(eulerangles)
% Mdc = MRIeuler2Mdc(eulerangles)
% 
% Returns the matrix of direction cosines given the Euler angles.
% Not sure if this works.
%
% $Id: MRIeuler2Mdc.m,v 1.1 2005/10/12 05:40:27 greve Exp $

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