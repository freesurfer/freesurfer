function Fx = binomialcdf(x,ntrials,theta)
% Fx = binomialcdf(x,ntrials,theta)
%
% Binomial cumulative probability distribution function. Returns the
% probability of seeing x or fewer events in ntrials when the expected
% value of x is ntrials*theta (ie, the likelihood of event X is
% theta).
%
% x must be a list integers where 0 < x < ntrials.
% 0 < theta < 1
%
% See also binomialpdf and binomialconf
%
% $Id: binomialcdf.m,v 1.1 2004/06/11 17:18:38 greve Exp $

Fx = [];

if(nargin ~= 3)
  fprintf('Fx = binomialcdf(x,ntrials,theta)\n');
  return;
end
nx = length(x);

% check for x all integers > 0 %
m = max(x-round(x));
if(m ~= 0)
  fprintf('ERROR: all x must be integers > 0\n');
  return;
end

% check for all x < n %
nover = length(find(x >= ntrials));
if(nover > 0)
  fprintf('ERROR: all x must be < ntrials\n');
  return;
end

% check 0 < theta < 1 %
if(theta <= 0 | theta >= 1)
  fprintf('ERROR: theta = %g, must be between 0 and 1\n',theta);
  return;
end

xmax = max(x);
xpdf = [1:xmax];
pdf = binomialpdf(xpdf,ntrials,theta);

% Use cumsum because distribution is discrete
cdf = cumsum(pdf);

Fx = cdf(x);

return;