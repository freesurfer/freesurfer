function [im, neigused] = invillcond(m,pvemin)
% [im neigused]= invillcond(m,pvemin)

im = [];
if(nargin ~= 1 & nargin ~= 2)
  msg = 'USAGE: [im neigused]= invillcond(m,pvemin)';
  qoe(msg);error(msg);
end

if(nargin == 1)
  pvemin = 100;
end

% Decompose %
[u s v] = svd(m);

% Vector of Eigenvalues %
meig = diag(s);
neig = length(meig);

% Cumulative Percent Variance Explained
pve = 100*cumsum(meig)/sum(meig);

% Pickout enough eigenvalues to explain the min desired %
neigused = length(find(pve < pvemin));

% regularize %
mineig = meig(neigused);
meig2 = meig;
meig2(neigused+1:neig) = mineig;

invs = diag(meig2.^(-1));

im = u*invs*v'; %'


return;
