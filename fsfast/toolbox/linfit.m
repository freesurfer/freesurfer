function [slope, offset] = linfit(y,x)
%
% [slope offset] = linfit(y,x)
%
% Solves the equation:
%   y = slope * x + offset
%
% $Id: linfit.m,v 1.1 2003/03/04 20:47:41 greve Exp $

if(nargin ~= 2)
  msg = 'USAGE: [slope offset] = linfit(y,x)';
  qoe(msg);error(msg);
end 

if(length(y) ~= length(x))
  msg = 'Length of x and y must be the same';
  qoe(msg);error(msg);
end

x = reshape1d(x);
y = reshape1d(y);

x2 = [x ones(size(x))];

m = pinv(x2)*y;

slope  = m(1);
offset = m(2);

return;