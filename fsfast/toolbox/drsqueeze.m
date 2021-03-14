function [y, xthresh, nclip, psqueeze] = drsqueeze(x,thr)
%
% y = drsqueeze(x)
%
% Squeezes dynamic range by clipping a fraction (thr) of the
% elements of greatest value.
%
%
%


%
% drsqueeze.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

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
