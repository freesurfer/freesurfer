function M1 = vox2ras_0to1(M0)
% M1 = vox2ras_0to1(M0)
%
% Converts a 0-based vox2ras matrix to 1-based, ie,
% Pxyz = M0*[c r s 1]' = M1*[(c+1) (r+1) (s+1) 1]'
%
% $Id: vox2ras_0to1.m,v 1.1 2002/12/13 23:45:31 greve Exp $


M1 = [];

if(nargin ~= 1)
  fprintf('M1 = vox2ras_0to1(M0)\n');
  return;
end

Q = zeros(4);
Q(1:3,4) = ones(3,1);

M1 = inv(inv(M0)+Q);

return;







