function acor = fast_cvm2acor(cvm) 
%
% acor = fast_cvm2acor(cvm) 
%
% Computes the autocorrelation function from the 
% covariance matrix

acor = [];

if(nargin ~= 1)
  fprintf('USAGE: acor = fast_cvm2acor(cvm)\n');
  return;
end

ntrs = size(cvm,1);

acor = zeros(ntrs,1);
for n = 1:ntrs, 
  acor(n) = mean(diag(cvm,n-1));
end

acor = acor/acor(1);

return;
