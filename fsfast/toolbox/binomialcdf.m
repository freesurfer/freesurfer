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
%


%
% binomialcdf.m
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
nover = length(find(x > ntrials));
if(nover > 0)
  fprintf('ERROR: all x must be <= ntrials\n');
  ntrials
  x
  return;
end

% check 0 < theta < 1 %
if(theta <= 0 | theta >= 1)
  fprintf('ERROR: theta = %g, must be between 0 and 1\n',theta);
  return;
end

xmax = max(x);
xpdf = [0:xmax];
pdf = binomialpdf(xpdf,ntrials,theta);

% Use cumsum because distribution is discrete
cdf = cumsum(pdf);

Fx = cdf(x+1);

% Can be > 1 because of round off
ind = find(Fx > 1);
Fx(ind) = 1;

return;
