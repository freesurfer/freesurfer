function ipeak = pkpulse(pulse,thresh)
% ipeak = pkpulse(pulse,thresh)
% 
% Returns the indices of the peaks in pulse.  Algorithm: convert pulse 
% into a vector of 0s and 1s by thresholding it at thresh.  If thresh is
% at the right level, this will create boxcar train of 1s and 0s.  This is 
% differentiated so that the begining of the boxcar is 1 and the end is -1.
% The peak is between these two limits is returned for each boxcar.
%
% A plot of just the peaks can be obtained by:
%   pulsepeaks = zeros(size(pulse));
%   pulsepeaks(ipeak) = pulse(ipeak);
%
% $Id: pkpulse.m,v 1.1 2003/03/04 20:47:41 greve Exp $

if(nargin ~= 2) 
  msg = 'USAGE: ipeak = pkpulse(pulse,thresh)';
  qoe(msg); error(msg);
end

% Threshold and differentiate
dth = diff(pulse > thresh);

% Get beginning (a) and ending (b) of boxcar %
a = find(dth ==  1);
b = find(dth == -1);
clear dth;

% Make sure there are the same number of beginnings and endings %
la = length(a);
lb = length(b);
nmin = min(la,lb);
if(la > nmin) a = a(1:nmin); end
if(lb > nmin) b = b(1:nmin); end

% Go through each boxcar and find the peak %
ipeak = zeros(1,nmin);
for n = 1:nmin
  rng = [a(n):b(n)];
  [m i] = max(pulse(rng));
  ipeak(n) = a(n) + i - 1;
end


return
