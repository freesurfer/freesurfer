function M = fast_blockdiag2(varargin)
% M = fast_blockdiag2(M1,M2,...)
%
% Assmbles matrix inputs into a block diagonal matrix.
% 
% See also fast_blockdiag
%
% $Id: fast_blockdiag2.m,v 1.1 2004/04/12 01:03:31 greve Exp $

M = [];
if(nargin == 0)
  fprintf('M = fast_blockdiag2(M1,M2,...)\n');
  return;
end

ncols = 0;
nrows = 0;
for n = 1:nargin
  Mn = varargin{n};
  nrows = nrows + size(Mn,1);
  ncols = ncols + size(Mn,2);
end

M = zeros(nrows,ncols);

r1 = 1;
c1 = 1;
for n = 1:nargin
  Mn = varargin{n};
  [nr nc] = size(Mn);
  r2 = r1 + nr - 1;
  c2 = c1 + nc - 1;
  M(r1:r2,c1:c2) = Mn;
  r1 = r2 + 1;
  c1 = c2 + 1;
end

return;







return;