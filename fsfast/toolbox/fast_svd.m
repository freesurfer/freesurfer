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
% 
% In any case, y = u*s*v';
%
% $Id: fast_svd.m,v 1.1 2004/08/18 20:05:40 greve Exp $
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
  u = y*(v*inv(s));
  uss2 = sqrt(sum(u.^2));
  u = u./repmat(uss2,[nr 1]);
  return;
end

% only gets here if(nr < nc)

M = y*y';
[u s blah] = svd(M);
s = sqrt(s);
v = y'*(u*inv(s));
vss2 = sqrt(sum(v.^2));
v = v./repmat(vss2,[nc 1]);

return;



