function [mse, MrExp] = fast_ar1w_mse(alpha,rho,Mr,R)
% [mse MrExp] = fast_ar1w_mse(alpha,rho,Mr,R)

mse = [];

if(nargin ~= 4)
  fprintf('mse = fast_ar1w_mse(alpha,rho,Mr,R)\n');
  return;
end

nf = size(R,1);
acf = fast_ar1w_acf(alpha,rho,nf);
MrExp = R*toeplitz(acf)*R;
E = Mr-MrExp;

%W = toeplitz(.99.^[0:nf-1]);
%E = W*E;

mse = sum(reshape1d(E).^2)/(nf*nf);

return;
