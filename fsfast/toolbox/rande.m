function r = rande(rdim,order)
%
% r = rande(dim,<order>)
%
% Generates matrix of given dimension whose elements 
% have an erlang distribution. The expected value
% is 1. The variance is 1/r. The stddev is 1/sqrt(r).
%
% If the order is unspecified, defalts to 1.
% An order of 1 is an exponential distribution.
%
% r = rande([5 50],3);
%
% See also erlang.m
%
% $Id: rande.m,v 1.3 2004/02/12 04:14:34 greve Exp $

if(exist('rdim')  ~= 1) rdim  = 1; end
if(exist('order') ~= 1) order = 1; end

r = 0;
for n = 1:order
  r = r + -log(rand(prod(rdim),1));
end
r = r/order; % make avg=1
r = reshape(r, [rdim 1]);

return;
