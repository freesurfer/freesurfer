function twarp = fast_timewarp(t,indpeak)
% twarp = fast_timewarp(t,indpeak)
%
% Change time so that the duration between peaks is constant and equal
% to the period=mean(diff(t(indpeak))). This is used for cases where
% there is a quasiperiodic process (like heat beat or respiration)
% where the time between peaks varies from cycle to cycle. This can be
% used in RETROICOR. For a truly periodic process, twarp=tresp (or
% very close to it). The time between peaks of the twarp should have
% the same/close mean of that for t but the stddev should be zero (or
% close to it). When the quasiperiodic process is plotted with
% twarp instead of t, then it should appear periodic.

%
% fast_timewarp.m
%
% Original Author: Douglas Greve
% Jan 29, 2019
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

npeaks = length(indpeak);
period = mean(diff(t(indpeak)));

nt = length(t);
dt = t(2)-t(1);
twarp = zeros(size(t));

% Handle before the first peak. This is just linear
k1 = 1;
k2 = indpeak(1)-1;
twarp(k1:k2) = t(k1:k2);

% Handle the time between the peaks.
for nthpeak = 2:npeaks
  k1 = indpeak(nthpeak-1);
  k2 = indpeak(nthpeak)-1;
  t1 = t(k1);
  t2 = t(k2+1);
  twarp(k1:k2) = period*(t(k1:k2)-t1)/(t2-t1) + (twarp(k1-1) + dt);
end

% Handle the time after the last peak
if(indpeak(end) < nt)
  k1 = indpeak(end);
  k2 = nt;
  t1 = t(k1);
  twarp(k1:k2) = twarp(k1-1) + [0:k2-k1]*dt;
end

return









