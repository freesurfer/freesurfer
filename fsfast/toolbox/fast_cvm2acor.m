function acor = fast_cvm2acor(cvm,biased) 
%
% acor = fast_cvm2acor(cvm) 
%
% Computes the autocorrelation function from the 
% covariance matrix

acor = [];

if(nargin ~= 1 & nargin ~= 2)
  fprintf('USAGE: acor = fast_cvm2acor(cvm,<biased>)\n');
  return;
end

if(exist('biased') ~= 1) biased = []; end
if(isempty(biased)) biased = 0; end

ntrs = size(cvm,1);

acor = zeros(ntrs,1);
for n = 1:ntrs, 
  acor(n) = mean(diag(cvm,n-1));
end

acor = acor/acor(1);

if(biased)
  acor = acor.*([ntrs:-1:1]');
  acor = acor/max(acor);
end


return;
