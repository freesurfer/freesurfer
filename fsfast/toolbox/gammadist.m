function y = gammadist(x,avg,stddev)
% Gamma Distribution
%
% y = gammadist(x,avg,stddev)
%
% This probably does not work.

if(nargin ~= 3) 
  msg = 'USAGE: y = gammadist(x,avg,stddev)'
  error(msg);
end

z = (x-avg)/stddev;
y = (z.^2) .* exp(-z);
ind = find(z<0);
y(ind) = 0;
y = y/sum(y);

return
