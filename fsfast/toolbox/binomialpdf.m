function fx = binomialpdf(x,n,theta)
% fx = binomialpdf(x,n,theta)
%
% Binomial probability distribution function. Returns the probability
% of seeing x events in n trials when the expected value of x is
% n*theta.
%
% x must be a list integers where 0 < x < n.
% 0 < theta < 1
%
% $Id: binomialpdf.m,v 1.1 2003/08/27 03:55:00 greve Exp $

fx = [];

if(nargin ~= 3)
  fprintf('fx = binomialpdf(x,n,theta)\n');
  return;
end
nx = length(x);

% check for x all integers > 0 %
m = max(x-round(x));
if(m ~= 0)
  fprintf('ERROR: all x must be integers > 0\n');
  return;
end

% check for all x < n %
nover = length(find(x >= n));
if(nover > 0)
  fprintf('ERROR: all x must be < n\n');
  return;
end

% check 0 < theta < 1 %
if(theta <= 0 | theta >= 1)
  fprintf('ERROR: theta = %g, must be between 0 and 1\n',theta);
  return;
end

% compute permuation for each x %
for nthx = 1:nx
  xx = x(nthx);
  %p(nthx) = factorial(n) / (factorial(xx) * factorial(n-xx) );
  p(nthx) = permutation(n,xx);
end

% finally compute pdf %
fx = p .* (theta.^x) .* ((1-theta).^(n-x));

return;

