function [pvs, u, ds] = fast_pvs(y)
% [pvs, u, s] = fast_pvs(y)
% 
% Computes percent variance spanned by each of the eigenvectors
%

pvs = [];
u = [];
ds = [];

if(nargin ~= 1)
  fprintf('[pvs, u, s] = fast_pvs(y)\n');
  return;
end

My = y*y'; %'
[u s v] = svd(My);
ds = diag(s);
pvs = 100*ds/sum(ds);

return;
