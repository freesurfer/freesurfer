function p = permutation(n,x)
% p = permutation(n,x);
% 
% Permutation of n and x:
% p = (n) =    n!
%     (x)   x!(n-x)!
%
% n is an integer > 0
% x is an integer > 0 and < n
%
% Uses logarithms to avoid overflow.
%
% $Id: permutation.m,v 1.1 2003/08/27 03:55:00 greve Exp $

p = [];
if(nargin ~= 2)
  fprintf('p = permutation(n,x)\n');
  return;
end

if(round(n) ~= n)
  fprintf('ERROR: n = %g, an integer\n',n);
  return;
end

if(round(x) ~= x)
  fprintf('ERROR: x = %g, an integer\n',x);
  return;
end

if(n <= 0)
  fprintf('ERROR: n = %d, must be > 0\n',n);
  return;
end

if(x <= 0 | x >= n)
  fprintf('ERROR: x = %d, must be > 0 and < n\n',x);
  return;
end

% Use logarithms to avoid overflows %
lp = sum(log(1:n)) - sum(log(1:x)) - sum(log(1:n-x));
p = exp(lp);

return;






