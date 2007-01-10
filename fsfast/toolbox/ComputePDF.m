function [pdf, xpdf, nxhist, cdf] = ComputePDF(x,xmin,xmax,xdelta)
%
% [pdf, xpdf, nxhist, cdf] = ComputePDF(x,xmin,xmax,xdelta)
%
% Creates a histogram nxhist with bins uniformly distributed 
% between xmin and xmax. Each bin is width xdelta. xpdf is
% the right edge of each bin. Using the right edge makes it
% possible to compute the CDF from cumsum(pdf).
%
% Ignores anything outside the range [xmin-xdelta xmax+xdelta]


%
% ComputePDF.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:29 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

pdf = [];
xpdf = [];
nxhist = [];
cdf = [];

if(nargin ~= 4)
  fprintf('[pdf, xpdf, nxhist] = ComputePDF(x,xmin,xmax,xdelta)\n');
  return;
end

% Get the number in the list before excluding out-of-range values
nlist = length(x);

% Ignore anything outside the range xmin-xdelta < x < xmax+xdelta %
ind = find(x >= (xmin-xdelta));
x = x(ind);
ind = find(x <= (xmax+xdelta));
x = x(ind);

nbins = round((xmax-xmin)/xdelta);

% Create a list of where the center of the bins are %
xcenter = xmin + ([0:nbins-1] + 0.5)*xdelta;
nxhist = hist(x,xcenter);
pdf = nxhist/nlist;

% Create a list of the right edge of the bin. This
% will work with cumsum().
xpdf = xmin+([1:nbins])*xdelta;

% Cumulative distribution
cdf = cumsum(pdf);

return;
