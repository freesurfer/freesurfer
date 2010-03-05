function AoB = fast_orthog(A,B)
% AoB = fast_orthog(A,B)
% Orthonalize A wrt B
% AoB = A - B*(inv(B'*B)*(B'*A));
%
% $Id: fast_orthog.m,v 1.1 2010/03/05 20:01:19 greve Exp $

AoB = [];
if(nargin ~= 2)
  fprintf('AoB = fast_orthog(A,B)\n');
  return;
end

AoB = A - B*(inv(B'*B)*(B'*A));

return;

