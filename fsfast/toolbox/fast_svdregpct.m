function [mreg, mdim] = fast_svdregpct(m,pct)
% [mreg, mdim] = fast_svdregpct(m,pct)
%
% Regularizes a matrix by choosing enough eigenvectors to span pct
% percent of the variance.
%
% $Id: fast_svdregpct.m,v 1.1 2003/09/18 02:24:22 greve Exp $

mreg = [];
if(nargin ~= 2)
  fprintf('[mreg, mdim] = fast_svdregpct(m,pct)\n');
  return;
end

[u s v] = svd(m);
ds = diag(s);
pvs = 100*ds/sum(ds);
cpvs = cumsum(pvs);

mdim = min(find(cpvs > pct));
ds2 = ds;
ds2(mdim:end) = ds2(mdim);

mreg = u * diag(ds2) * v';

return;

















