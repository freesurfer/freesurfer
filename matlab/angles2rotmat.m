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
% $Id: angles2rotmat.m,v 1.4 2011/03/02 00:04:12 nicks Exp $

%
% angles2rotmat.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.4 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
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



