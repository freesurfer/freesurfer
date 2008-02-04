function [g2 G2] = mri_kurtosis(x,dim)
% g2 = mri_kurtosis(x,<dim>)
% Biased estimator of kurtosis.

g2 = [];
if(nargin ~= 1 & nargin ~= 2)
  fprintf('g2 = mri_kurtosis(x,<dim>)\n');
  return;
end

if(nargin == 1) dim = []; end
if(isempty(dim)) dim = 1; end

szx = size(x);
n = szx(dim);

xmean = mean(x,dim);
tmp = ones(1,length(szx));
tmp(dim) = n;
xmeanRep = repmat(xmean,tmp);

m4 = n*sum((x - xmeanRep).^4,dim);
m2 = sum((x - xmeanRep).^2,dim);
%g2 = m4./(m2.^2) - 3;
g2 = m4./(m2.^2);

b1 = (n+1)*(n-1)/((n-2)*(n-3));
b2 = ((n-1).^2)/((n-2)*(n-3));
%G2 = b1*(m4./(m2.^2)) - 3*b2;
G2 = b1*(m4./(m2.^2));

return;
