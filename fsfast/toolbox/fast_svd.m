function [u,s,v] = fast_svd(y)
% [u s v] = fast_svd(y)
% 
% Computes efficient SVD when the number of rows and columns are
% not the same. It is efficient in the sense that only the minimum
% number of eigen components are computed.
%
% u will have dimension nrows by nmin
% v will have dimension ncols by nmin
% s will have dimension nmin by nmin
%     where nmin = min(nrows,ncols)
%     nmin can also be determined as the number of non-zero singvals
% 
% In any case, y = u*s*v';
%
% $Id: fast_svd.m,v 1.2 2004/08/19 00:53:07 greve Exp $
%

u=[];
s=[];
v=[];

if(nargin ~= 1)
  fprintf('[u s v] = fast_svd(y)\n');
  return;
end

[nr nc] = size(y);

if(nr == nc) 
  [u s v] = svd(y); 
  return;
end

if(nr > nc)
  M = y'*y;
  [v s blah] = svd(M);
  s = sqrt(s);
  ds = diag(s);
  ns = length(find(ds > ds(1)/1e6));
  nn = 1:ns;
  s = s(nn,nn);
  v = v(:,nn);
  u = y*(v*inv(s));
  uss2 = sqrt(sum(u.^2));
  u = u./repmat(uss2,[nr 1]);
  return;
end

% only gets here if(nr < nc)

M = y*y';
[u s blah] = svd(M);
s = sqrt(s);
ds = diag(s);
ns = length(find(ds > ds(1)/1e6));
nn = 1:ns;
s = s(nn,nn);
u = u(:,nn);
v = y'*(u*inv(s));
vss2 = sqrt(sum(v.^2));
v = v./repmat(vss2,[nc 1]);

return;



