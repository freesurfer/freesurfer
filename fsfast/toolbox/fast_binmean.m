function [m, s, n, xmin, xmax] = fast_binmean(x,bincenter,y)
% [m, s, n, xmin, xmax] = fast_binmean(x,bincenter,<y>)
%
% m = mean, s = stddev, n = number in each bin
% If there are no data points in a bin, the mean and stddev
% are set to zero.
%
% If y is present, then it must be the same length as x, and
% the return values are that of the mean and stddev of y
% within the bins instead of x.
%
% xmin and xmax are the boundaries of the bins.
%

m = [];
s = [];

if(nargin ~= 2 & nargin ~= 3)
  fprintf('USAGE: [m, s] = fast_binmean(x,bincenter,<y>)\n');
  return;
end

if(nargin ~= 3) y = x; end;
if(length(x) ~= length(y))
  fprintf('ERROR: x and y have different lengths\n');
  return;
end

dbin = bincenter(2) - bincenter(1);

nbins = length(bincenter);
xmin = bincenter - dbin/2;
xmax = bincenter + dbin/2;

for bin = 1:nbins
  if(bin ~= nbins)
    ind = find(x >= xmin(bin) & x < xmax(bin));
  else
    ind = find(x >= xmin(bin) & x <= (xmax(bin)+10^(-10)) );
  end
  n(bin) = length(ind);
  if(~isempty(ind))
    m(bin) = mean(y(ind));
    s(bin) = std(y(ind));
  else
    m(bin) = 0;
    s(bin) = 0;
  end
end

return;
