function [y, xthresh, nclip, psqueeze] = drsqueeze(x,thr)
%
% y = drsqueeze(x)
%
% Squeezes dynamic range by clipping a fraction (thr) of the
% elements of greatest value.
%
%
% $Id: drsqueeze.m,v 1.1 2003/03/04 20:47:34 greve Exp $

xsize = size(x);

x = reshape1d(x);
nx = length(x);
nb = floor(nx/20);
if(nb > 200) nb = 200; end

xmin = min(x);
xmax = max(x);

[xhist xbin] = hist(x,nb);

cdf = cumsum((xhist/nx));

ithresh = min(find(cdf>thr));
xthresh = xbin(ithresh);

isupthresh = find(x>xthresh);
y = x;
y(isupthresh) = xthresh;
y = reshape(y,xsize);

nclip = length(isupthresh);
psqueeze = (xthresh-xmin)/(xmax-xmin);

return;
