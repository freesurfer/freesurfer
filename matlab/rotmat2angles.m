function [angles ambflag]= rotmat2angles(R)
% [angles ambflag]= rotmat2angles(R)
%
% Converts rotation matrix to euler angles (in radians)
% angles is a 3x1 vector in radians
% angles(1) - pitch - rotation about x or LR (gamma)
% angles(2) - yaw   - rotation about y or AP (beta)
% angles(3) - roll  - rotation about z or SI (alpha)
%
% ambflag=1 indicates some ambiguity (eb, cos(beta)=1)
% There may also be amibuity if the rotation matrix
% was created with angles > pi/2.
%
% See also: angles2rotmat
% Ref: Craig, Intro to Robotics
%
% $Id: rotmat2angles.m,v 1.2 2007/07/25 19:30:01 greve Exp $

%
% rotmat2angles.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/07/25 19:30:01 $
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

angles = [];
ambflag = [];
if(nargin ~= 1)
  fprintf('angles = rotmat2angles(R)\n');
  return;
end

r31 = R(3,1);
beta = asin(-r31);

% This threshold is surprisingly large
if(abs(1-abs(r31)) > 1e-4)
  alpha = atan2(R(2,1),R(1,1));
  gamma = atan2(R(3,2),R(3,3));
  ambflag = 0;
else
  % beta = +/-pi/2 = +/-90deg -- ambiguous?
  alpha = 0;
  gamma = atan2(-R(1,2),R(2,2));
  ambflag = 1;
  %printf('here\n');
end

angles = [gamma beta alpha ]';

return;








