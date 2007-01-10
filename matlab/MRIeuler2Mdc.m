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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
%    $Revision: 1.2 $
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
