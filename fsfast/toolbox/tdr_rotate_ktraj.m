function [kxr, kyr] = tdr_rotate_ktraj(kx,ky,theta)
% [kxr, kyr] = tdr_rotate_ktraj(kx,ky,theta)
%
% Rotate a k-space trajectory counter-clockwise by
% by theta radians about the center of k-space. This 
% can be used to simulate a propeller sequence.
%
% $Id: tdr_rotate_ktraj.m,v 1.1 2004/01/23 20:17:11 greve Exp $

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









