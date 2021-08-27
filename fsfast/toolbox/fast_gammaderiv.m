function [dh, t] = fast_gammaderiv(delta, tau, TER, TPreStim, TimeWindow)
%
% [dh t] = fast_gammaderiv(delta, tau, TER, TPreStim, TimeWindow)
%
% Computes the derivative of the gamma function (from Dale and 
% Buckner, 1997):
%
% h(t>delta)  = ((t-delta)/tau)^2 * exp(-(t-delta)/tau)
% h(t<=delta) = 0;
%
% The waveform is computed by sampling the continuous derivative (ie,
% it is not computed by computing the gamma function and then doing
% a diff). The waveform is scaled so that the peak = 1.0
%
% Sample parameters (from D&B): delta = 2.25 sec, tau = 1.25 sec
%
% See also fast_gamma.
%
%
%


%
% fast_gammaderiv.m
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

dh = [];

if(nargin ~= 5)
  msg = 'USAGE: [h t] = fast_gamma(delta, tau, TER, TPreStim, TimeWindow)\n';
  qoe(msg);error(msg);
end

NWindow = round(TimeWindow/TER);
t = TER * [0:NWindow-1] - TPreStim;
t = reshape1d(t); % Column vector

r = (t - delta)./tau ;
dh = ((2-r) .* r .* exp(-r))/tau;
i0 = find(t<delta);
dh(i0) = zeros(size(t(i0)));

% scale max to 1 %
dh = dh/max(dh);

return;
