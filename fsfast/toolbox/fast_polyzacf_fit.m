function [coeff, acffit, X] = fast_polyzacf_fit(acf,order,maxlag,wflag)
% [coeff, acffit, X] = fast_polyzacf_fit(acf,order,maxlag,wflag)
%
% Fits the autocorrelation function acf to a polynomial that
% forces acffit to be zero at maxlag.
%
% maxlag only the first maxlag components of acf are fit.
%
%  acffit = X*coeff;
%  acffit = [ones(1,size(acf,2)); acffit];
%
%


%
% fast_polyzacf_fit.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
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

coeff = [];
acffit = [];
X = [];

if(nargin < 2 | nargin > 4)
  msg = 'USAGE [coeff, acffit] = fast_polyzacf_fit(acf,order,<maxlag>,<wflag>)';
  qoe(msg); error(msg);
end

nacf = size(acf,1);
nv = size(acf,2);
if(exist('maxlag') ~= 1) maxlag = nacf; end
if(isempty(maxlag)) maxlag = nacf; end

if(exist('wflag') ~= 1) wflag = 0; end
if(isempty(wflag)) wflag = 0; end

% remove zero lag everything beyond maxlag %
acf = acf(2:maxlag,:);

rhomin = 0;
rhomax = 0.9;
drho = (rhomax-rhomin)/(order-1);
tau = [1:maxlag-1]';

X = [];
for n=1:order
  x = (max(tau)-tau).^n;
  x = x/std(x);
  X = [X x];
end
[u s v] = svd(X);
ds = diag(s);
cpvs = cumsum(ds)/sum(ds);
ind = find(cpvs<.99);
X = u(:,ind);

if(wflag ~= 0)
  X0 = X;
  w = ([maxlag-1:-1:1]).^wflag;
  w = w/max(w);
  W = diag(w);
  X   = W*X;
  acf = W*acf;
end

coeff = (inv(X'*X)*X')*acf;

if(nargout > 1)
  if(wflag ~= 0) X = X0; end
  acffit = X*coeff;
  acffit = [ones(1,nv); acffit];
  if(maxlag ~= nacf)
    acffit = [acffit; zeros(nacf-maxlag,nv)];
  end
end

return;






