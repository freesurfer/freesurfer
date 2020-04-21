function AoB = fast_orthog(A,B)
% AoB = fast_orthog(A,B)
% Orthonalize A wrt B
% AoB = A - B*(inv(B'*B)*(B'*A));
%

AoB = [];
if(nargin ~= 2)
  fprintf('AoB = fast_orthog(A,B)\n');
  return;
end

AoB = A - B*(inv(B'*B)*(B'*A));

return;

