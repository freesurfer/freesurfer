function y = fast_monotonize(x)
% y = fast_monotonize(x)
%
% Makes columns of x monotonically decreasing
%
% $Id: fast_monotonize.m,v 1.1 2004/04/11 07:24:17 greve Exp $

if(nargin ~= 1)
  fprintf('y = fast_monotonize(x)\n');
  return;
end

nf = size(x,1);

y = x;
for n = 2:nf
  ind = find(y(n,:) > y(n-1,:));
  y(n,ind) = y(n-1,ind) ;
end

return;