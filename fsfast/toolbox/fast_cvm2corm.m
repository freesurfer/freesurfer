function cormtx = fast_cvm2corm(cvm)
% cormtx = fast_cvm2corm(cvm)
% Converts a covariance matrix into a correlation coefficient
% matrix. It computes the correlation coefficient for compoent
% i,j as cvm(i,j) / sqrt(cvm(i,i) * cvm(j,j))

cormtx = [];

if(nargin ~= 1)
  msg = 'USAGE: cormtx = fast_cvm2corm(cvm)';
  fprintf('%s\n',msg);
  return;
end

d = diag(cvm); % n X 1
n = length(d);

m = sqrt(repmat(d, [1 n]) .* repmat(d', [n 1])); %'

cormtx = cvm ./ m;

return;
