function [m, i] = mar(x)
%
% [m i] = mar(x)
%
% [m i] = max(abs(reshape1d(x)))
%

if(nargin ~= 1)
  msg = 'USAGE: [m i] = mar(x)';
  qoe(msg);error(msg);
end

n = prod(size(x));

[m i] = max(abs(reshape(x,[n 1])));

return;
