function [edge, bincenter, binmap] = fast_histeq(y,nbins,r)
% [edge, bincenter, binmap] = fast_histeq(y,nbins,<r>)
%
% computes the bin edges that will result in an equal number of
% samples of y in each bin. The first step in the algorithm is to
% compute a histogram of y; the number of bins is nbins*r. If r is
% unspecfied, it defaults to 20. The larger r is, the better the final
% result, but the more data needed.
%
% The size of edge will be nbins+1, where the first bin will be
% between edge(1) and edge(2), etc.
%
%
% To check:
% y = randn(10000,1);
% edge = fast_histeq(y,100);
% nk = histc(y,edge);
% plot(nk(1:end-1));
% plot should have approx 1000 = 10000/nbins at each entry
%
% $Id: fast_histeq.m,v 1.3 2004/06/11 17:21:24 greve Exp $

edge = [];

if(nargin ~= 2 & nargin ~= 3)
  fprintf('edge = fast_histeq(y,nbins,<r>)\n');
  return;
end

if(exist('r') ~= 1) r = []; end
if(isempty(r)) r = 20; end

ny = length(y);
nbins0 = r*nbins;

if(ny < nbins0)
  fprintf('ERROR: not enough data points for %d bins (r=%d)\n',nbins,r);
  return;
end

[nperbin bincenter] = hist(y,nbins0);
p = nperbin/ny;
cp = cumsum(p);
binwidth = mean(diff(bincenter));

edge = zeros(nbins+1,1);
edge(1)   = bincenter(1) - binwidth/2;
edge(end) = bincenter(end) + binwidth/2;
for n = 1:nbins-1
  pedge = n/nbins;
  [m i] = min(abs(cp-pedge));
  edge(n+1) = bincenter(i) + binwidth/2;
end

if(nargout < 2) return; end
bincenter = (edge(1:end-1)+edge(2:end))/2;

if(nargout < 3) return; end

binmap = zeros(size(y));
for n = 1:nbins
  if(n == 1)
    ind = find(y <= edge(n+1) );
  elseif(n==nbins)
    ind = find(edge(n) < y);
  else
    ind = find(edge(n) < y & y <= edge(n+1) );
  end
  binmap(ind) = n;
end



%nk = histc(y,edge);
%plot(nk(1:end-1));
%keyboard

return;























