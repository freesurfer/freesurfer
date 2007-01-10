function [kxr, kyr] = tdr_rotate_ktraj(kx,ky,theta)
% [kxr, kyr] = tdr_rotate_ktraj(kx,ky,theta)
%
% Rotate a k-space trajectory counter-clockwise by
% by theta radians about the center of k-space. This 
% can be used to simulate a propeller sequence.
%
%


%
% tdr_rotate_ktraj.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:35 $
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

kxr = [];
kyr = [];

if(nargin ~= 3)
  fprintf('[kxr, kyr] = tdr_rotate_ktraj(kx,ky,theta)\n');
  return;
end

nx = prod(size(kx));
ny = prod(size(ky));
if(nx ~= ny)
  fprintf('ERROR: kx and ky have different lengths\n');
  return;
end

r = sqrt(kx.^2 + ky.^2);
phi = atan2(ky,kx);

kxr = r.*cos(phi+theta);
kyr = r.*sin(phi+theta);


return;









