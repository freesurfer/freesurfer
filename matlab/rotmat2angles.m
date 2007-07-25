function angles = rotmat2angles(R)
% angles = rotmat2angles(R)
%
% WARNING: this will fail when cos(angle(1)) = pi/2
%
% Converts rotation matrix to angles (in radians)
% angles(1) - pitch - rotation about x or LR
% angles(2) - yaw   - rotation about y or AP
% angles(3) - roll  - rotation about z or SI
%
% See also: angles2rotmat
%
% $Id: rotmat2angles.m,v 1.1 2007/07/25 06:23:26 greve Exp $

%
% rotmat2angles.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/07/25 06:23:26 $
%    $Revision: 1.1 $
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

angles = [];
if(nargin ~= 1)
  fprintf('angles = rotmat2angles(R)\n');
  return;
end

r32 = R(3,2);
alpha = -asin(r32);

r31 = R(3,1);
beta = asin(r31/cos(alpha));

r22 = R(2,2);
gamma = acos(r22/cos(alpha));


angles = [alpha beta gamma];

return;








