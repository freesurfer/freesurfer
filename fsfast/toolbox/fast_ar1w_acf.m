function acf = fast_ar1w_acf(alpha, rho, nf)
% acf = fast_ar1w_acf(alpha, rho, nf)
%
% Autocorrelation function of AR1+White Noise
%   acf(0) = 1
%   acf(n) = (1-alpha)*rho^n

acf = [];

if(nargin ~= 3)
  fprintf('acf = fast_ar1w_acf(alpha, rho, nf)\n');
  return;quit;
end

nn = [0:nf-1]'; %'
acf = (1-alpha)*rho.^nn;
acf(1) = 1;

return;
