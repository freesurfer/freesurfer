function M0 = vox2ras_1to0(M1)
% M0 = vox2ras_1to0(M1)
%
% Converts a 1-based vox2ras matrix to 0-based, ie,
% Pxyz = M0*[c r s 1]' = M1*[(c+1) (r+1) (s+1) 1]'
%
% See also: vox2ras_0to1
%
% $Id: vox2ras_1to0.m,v 1.1 2003/01/09 23:48:26 greve Exp $

M0 = [];

if(nargin ~= 1)
  fprintf('M0 = vox2ras_1to0(M1)');
  return;
end

Q = zeros(4);
Q(1:3,4) = -ones(3,1);

M0 = inv(inv(M1)+Q);

return;







