function M = fast_ar1mtx(rho,len)
% M = fast_ar1mtx(rho,len)
%
% Produces a toeplitz AR1 matrix
%

M = [];

if(nargin ~= 2)
  fprintf('USAGE: M = fast_ar1mtx(rho,len)\n');
  return;
end

M = [];
for c = 0:len-1;
  ind = abs([-c : len-c-1]);
  r = rho .^ ind;
  M = [M; r];
end

return;
