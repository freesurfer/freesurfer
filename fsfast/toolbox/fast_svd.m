function [u,s,v] = fast_svd(y,M)
% [u s v] = fast_svd(y,<M>)
% 
% Computes efficient SVD when the number of rows and columns are
% not the same. It is efficient in the sense that only the minimum
% number of eigen components are computed. If M is supplied, then
% must be M = y*y' or M = y'*y depending upon the size of y.
%
% u will have dimension nrows by nmin
% v will have dimension ncols by nmin
% s will have dimension nmin by nmin
%     where nmin = min(nrows,ncols)
%     nmin can also be determined as the number of non-zero singvals
% 
% In any case, y = u*s*v';
%
% $Id: fast_svd.m,v 1.3 2004/09/23 18:25:39 greve Exp $
%

u=[];
s=[];
v=[];

if(nargin < 1 | nargin > 2)
  fprintf('[u s v] = fast_svd(y,<M>)\n');
  return;
end

[nr nc] = size(y);

if(nr == nc) 
  [u s v] = svd(y); 
  return;
end

if(nr > nc)
  if(~exist('M','var')) M = y'*y; end
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

if(~exist('M','var')) M = y*y'; end
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



