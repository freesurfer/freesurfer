function y = erlang(x,avg,r)
% y = erlang(x,avg,r)
%
% Generates theoretical erlang distribution with average avg and order
% r. If avg and r are not specfied, both are assumed to be 1. The
% variance will be (avg^2)/r. The stddev = avg/sqrt(r).
%
% Note: it is not possible to control the average and variance
% separately.
%
% To test emperically:
%   y = rande(10000,5);
%   [h x] = hist(y,100);  h = h/max(h);
%   hest = erlang(x,1,5); hest = hest/max(hest);
%   plot(x,h,x,hest);
%
% mu = 1/avg;
% y = r*mu*((r*mu*x)^(r-1)) * exp(-r*mu*x) / (r-1)!
%
% See also rande.
%
% $Id: erlang.m,v 1.1 2004/02/12 04:11:30 greve Exp $

y = [];

if(nargin < 1 | nargin > 3)
  fprintf('y = erlang(x,avg,r)\n');
  return;
end

if(exist('avg')~=1) avg = 1; end
if(exist('r')~=1)  r = 1; end

y = zeros(size(x));
indgez = find(x >= 0);

mu = 1/avg;
rmu = r*mu;
y(indgez) = rmu*((rmu*x(indgez)).^(r-1)) .* exp(-rmu*x(indgez)) / factorial(r-1);

return;
  
  
  
  