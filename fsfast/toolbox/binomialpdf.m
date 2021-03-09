function fx = binomialpdf(x,ntrials,theta)
% fx = binomialpdf(x,ntrials,theta)
%
% Binomial probability distribution function. Returns the probability
% of seeing x events in ntrials when the expected value of x is
% ntrials*theta. (ie, the likelihood of event X is theta).
%
% x must be a list integers where 0 < x < ntrials.
% 0 < theta < 1
%
% See also binomialcdf and binomialconf
%
%


%
% binomialpdf.m
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

fx = [];

if(nargin ~= 3)
  fprintf('fx = binomialpdf(x,ntrials,theta)\n');
  return;
end
nx = length(x);

% check for x all integers > 0 %
m = max(x-round(x));
if(m ~= 0)
  fprintf('ERROR: all x must be integers > 0\n');
  return;
end

% check for all x <= ntrials %
nover = length(find(x > ntrials));
if(nover > 0)
  fprintf('ERROR: all x must be <= ntrials\n');
  return;
end

% check 0 < theta < 1 %
if(theta <= 0 | theta >= 1)
  fprintf('ERROR: theta = %g, must be between 0 and 1\n',theta);
  return;
end

% compute log of permuation for each x %
for nthx = 1:nx
  xx = x(nthx);
  %p(nthx) = factorial(n) / (factorial(xx) * factorial(n-xx) );
  %p(nthx) = permutation(ntrials,xx)
  lp(nthx) = permutation(ntrials,xx,1);
end

% finally compute pdf %
%fx = p .* (theta.^x) .* ((1-theta).^(ntrials-x));
lfx = lp + x.*log(theta) + (ntrials-x)*log(1-theta);
fx = exp(lfx);

return;

