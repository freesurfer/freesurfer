function [xlow, xhi, xmean, xstd] = binomialconf(n,theta,pct)
% [xlow, xhi, xmean, xstd] = binomialconf(ntrials,theta,pct)
%
% Computes the theoretical confidence intervals for a binomial
% distribution. For a given experiment with n trials where theta is
% the probability that the result of a trial is of type X, and x is
% the actual number of type X, then the expected value of x is xmean =
% ntrials*theta, the standard deviation of x is xstd =
% sqrt(ntrials*theta*(1-theta)), and pct percent of the time x will
% fall between xlow and xhi (inclusive). If p = x/ntrials, then pct
% percent of the time, p will fall between xlow/ntrials and
% xhigh/ntrials (inclusive). This range of p will not be symetrical
% around x/ntrials.
%
% pct is a percent (0-100).
% ntrials must be an integer greater than 0.
% 0 < theta < 1
% 0 < pct < 100.
%
% See also binomialpdf and binomialcdf.
%
%


%
% binomialconf.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:29 $
%    $Revision: 1.3 $
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

xmean = [];
xstd = [];
xlow = [];
xhi = [];

if(nargin ~= 3)
  fprintf('[xlow, xhi, xmean, xstd] = binomialconf(ntrials,theta,pct)\n');
  return;
end

if(pct <= 0 | pct >= 100)
  fprintf('ERROR: pct = %g, must be 0 < pct < 100\n',pct);
  return;
end

xmean = round(n*theta);
xstd  = sqrt(n*theta*(1-theta));

% Create a list of x that bracket the mean %
x = round([xmean-4*xstd:xmean+4*xstd]);
indok = find(x > 0 & x < n);
x = x(indok);

% Generate the pdf at each of those points.
fx = binomialpdf(x,n,theta);
if(isempty(fx)) return; end

% Compute the cdf
cdfx = cumsum(fx);

% Find where the cdf passes through the low and hi points
d = (100-pct)/2; % tail to be excluded on either side
p_low = d/100;
p_hi  = 1 - p_low;
[m indlow] = min(abs(cdfx-p_low));
[m indhi]  = min(abs(cdfx-p_hi));

% Record those values of x
xlow = x(indlow);
xhi  = x(indhi);

return;


