function h = fmri_hemodyn(t, delta, tau)
%
% function h = fmri_hemodyn(t, delta, tau)
%
% Simulates the hemodynaimc response using model
% from Dale and Buckner, 1997:
%
% h(t>delta)  = ((t-delta)/tau)^2 * exp(-(t-delta)/tau)
% h(t<=delta) = 0;
%
% The HDIR is scaled so that the peak = 1.0
%
% Sample parameters: delta = 2.25 sec, tau = 1.25 sec
%
% $Id: fmri_hemodyn.m,v 1.1 2003/03/04 20:47:39 greve Exp $
%

if(nargin ~= 3)
  msg = 'USAGE: h = fmri_hemodyn(t, delta, tau)';
  qoe(msg);error(msg);
end

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
h = ( ( r.^2) .* exp(-r) );
i0 = find(t<delta);
h(i0) = zeros(size(t(i0)));

% scale max to 1 %
h = h*exp(2.0)/4.0;

return;
