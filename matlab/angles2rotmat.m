function [R Rx Ry Rz] = angles2rotmat(angles)
% [R Rx Ry Rz] = angles2rotmat(angles)
%
% Convert 3 euler angles into a rotation matrix
%
% angles is a 3x1 vector in radians
% angles(1) - pitch - rotation about x or LR (gamma)
% angles(2) - yaw   - rotation about y or AP (beta)
% angles(3) - roll  - rotation about z or SI (alpha)
% R = Rz*Ry*Rx;
%
% See also: rotmat2angles
% Ref: Craig, Intro to Robotics
%
% $Id: angles2rotmat.m,v 1.3 2007/07/25 19:30:01 greve Exp $

%
% angles2rotmat.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/07/25 19:30:01 $
%    $Revision: 1.3 $
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

R  = [];
Rx = [];
Ry = [];
Rz = [];
if(nargin ~= 1)
  fprintf('R = angles2rotmat(angles)\n');
  return;
end

gamma = angles(1);
beta  = angles(2);
alpha = angles(3);

Rx = zeros(3,3);
Rx(1,1) = +1;
Rx(2,2) = +cos(gamma);
Rx(2,3) = -sin(gamma);
Rx(3,2) = +sin(gamma);
Rx(3,3) = +cos(gamma);

Ry = zeros(3,3);
Ry(1,1) = +cos(beta);
Ry(1,3) = +sin(beta);
Ry(2,2) = +1;
Ry(3,1) = -sin(beta);
Ry(3,3) = +cos(beta);

Rz = zeros(3,3);
Rz(1,1) = +cos(alpha);
Rz(1,2) = -sin(alpha);
Rz(2,1) = +sin(alpha);
Rz(2,2) = +cos(alpha);
Rz(3,3) = +1;

R = Rz*Ry*Rx;

return;



