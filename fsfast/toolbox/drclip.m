function [y, xthresh, nclip, psqueeze] = drclip(x,thr,clipval)
%
% [y, xthresh, nclip, psqueeze] = drclip(x,thr,<clipval>)
%
% Squeezes dynamic range by clipping a fraction (thr) of the
% elements of greatest value to clipval.
%
%
% $Id: drclip.m,v 1.1 2003/03/04 20:47:33 greve Exp $

if(nargin ~= 2 & nargin ~= 3)
  msg = 'USAGE: [y, xthresh, nclip, psqueeze] = drclip(x,thr,<clipval>)';
  qoe(msg); error(msg);
end

if(nargin == 2) clipval = 0; end

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
y(isupthresh) = clipval;
y = reshape(y,xsize);

nclip = length(isupthresh);
psqueeze = (xthresh-xmin)/(xmax-xmin);

return;
