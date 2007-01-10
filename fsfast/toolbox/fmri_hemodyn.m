function h = fmri_hemodyn(t, delta, tau, alpha)
%
% function h = fmri_hemodyn(t, delta, tau, <alpha>)
%
% Simulates the hemodynaimc response using model
% from Dale and Buckner, 1997:
%
% h(t>delta)  = ((t-delta)/tau)^alpha * exp(-(t-delta)/tau)
% h(t<=delta) = 0;
%
% The HDIR is scaled so that the continuous-time peak = 1.0,
% though the peak of the sampled waveform may not be 1.0.
%
% Sample parameters: delta = 2.25 sec, tau = 1.25 sec, alpha=2
%
% Be default, alpha = 2 (good for BOLD).
%
%
%


%
% fmri_hemodyn.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
%    $Revision: 1.4 $
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

if(nargin ~= 3 & nargin ~= 4)
  msg = 'USAGE: h = fmri_hemodyn(t, delta, tau, <alpha>)';
  qoe(msg);error(msg);
end

if(nargin == 3) alpha = 2; end

if(length(delta) ~= length(tau))
  msg = 'delta and tau dimensions are inconsistent';
  qoe(msg);error(msg);
end

nh = length(delta);
nt = length(t);

t     = reshape(t,     [nt 1]);
delta = reshape(delta, [1 nh]);
tau   = reshape(tau,   [1 nh]);

t     = repmat(t,     [1 nh]);
delta = repmat(delta, [nt 1]);
tau   = repmat(tau,   [nt 1]);

r = (t - delta)./tau ;
h = ( ( r.^alpha) .* exp(-r) );
i0 = find(t<delta);
h(i0) = zeros(size(t(i0)));

% Scale so that max of continuous function is 1.
% Peak will always be at (alpha.^alpha)*exp(-alpha)
peak = (alpha.^alpha)*exp(-alpha);
h = h/peak;

return;
