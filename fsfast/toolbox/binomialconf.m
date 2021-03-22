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

ntheta = length(theta);
if(ntheta > 1)
  fprintf('binomialconf: looping %d\n',ntheta);
  xlow  = zeros(ntheta,1);
  xhi   = zeros(ntheta,1);
  xmean = zeros(ntheta,1);
  xstd  = zeros(ntheta,1);
  for nththeta = 1:ntheta
    [xlow(nththeta) xhi(nththeta) xmean(nththeta) xstd(nththeta)] = ...
	binomialconf(n,abs(theta(nththeta)),pct);
  end
  return;
end

if(theta <= 0)
  xlow = 0;
  xhi  = n;
  xmean = 0;
  xstd = 100;
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


