function [R Rx Ry Rz] = angles2rotmat(angles)
% [R Rx Ry Rz] = angles2rotmat(angles)
%
% Convert 3 angles into a rotation matrix
%
% angles is 3x1 in radians
% angles(1) - rotation about x
% angles(2) - rotation about y
% angles(3) - rotation about z
% R = Rz*Ry*Rx; All 3x3 matrices
%
% $Id: angles2rotmat.m,v 1.1 2007/07/25 05:56:13 greve Exp $

%
% angles2rotmat.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/07/25 05:56:13 $
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

R  = [];
Rx = [];
Ry = [];
Rz = [];
if(nargin ~= 1)
  fprintf('R = angles2rotmat(angles)\n');
  return;
end

alpha = angles(1);
beta  = angles(2);
gamma = angles(3);

Rx = zeros(3,3);
Rx(1,1) = +1;
Rx(2,2) = +cos(alpha);
Rx(2,3) = +sin(alpha);
Rx(3,2) = -sin(alpha);
Rx(3,3) = +cos(alpha);

Ry = zeros(3,3);
Ry(1,1) = +cos(beta);
Ry(1,3) = -sin(beta);
Ry(2,2) = +1;
Ry(3,1) = +sin(beta);
Ry(3,3) = +cos(beta);

Rz = zeros(3,3);
Rz(1,1) = +cos(gamma);
Rz(1,2) = +sin(gamma);
Rz(2,1) = -sin(gamma);
Rz(2,2) = +cos(gamma);
Rz(3,3) = +1;

R = Rz*Ry*Rx;

return;



