function [X, a, Xfir] = fast_sched2Xgamma(tPres,ntrs,TR,delay,dispersion,bcw,W)
%
% [X, a] = fast_sched2Xgamma(tPres,ntrs,TR,delay,dispersion,bcw,W)
%
% Creates a design matrix (aka stimulus convolution matrix) modeling
% the hemodynamic response as a single gamma function. The size of 
% the matrix X will be ntrs-by-1. a is the assumed form.
%
% tPres - list of presentation times (ie, the schedule) in seconds 
% for one event type. Time=0 is defined as the time the first stored 
% image was collected (ie, ignore discarded acquisitions). An event
% before t=0 will have an effect on the matrix if its Post Stimulus 
% Window encompasses t=0; otherwise it is ignored. Events found
% after the termination of data collection are ignored (a warning 
% is printed). If tPres is empty, a matrix of zeros is returned. 
%
% ntrs - number of functional volumes collected (do not include 
% discarded acquisitions or prescan time).
%
% TR - TR of the exmperiment (ie, time between acquisitions of 
% functional volumes) in seconds.
%
% delay - delay in the gamma function
% dispersion - time constant
% bcw - convolve with a boxcar of width bcw (seconds)
%
% W - presentation weighting. Length should be NPresentations. The
% matrix entry for each presentation is given the value W(n) instead
% of 1. Ignored if W=[].
%
% See also: fast_gamma, fast_sched2Xfir, fast_sched2Xgammaderiv
%
%


%
% fast_sched2Xgamma.m
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

X = [];
a = [];

if(nargin ~= 7)
  msg = 'USAGE: [X, a] = fast_sched2Xgamma(tPres,ntrs,TR,delay,dispersion,bcw,W)';
  fprintf('%s\n',msg);
  return;
end

if(isempty(bcw)) bcw = 0; end

% Compute an "appropriate" time window %
TimeWindow = delay + 10*dispersion;

% Compute an "appropriate" temporal resolution
Tres = TR/10;

% Recompute time window as mult of Tres %
TimeWindow = Tres*round(TimeWindow/Tres);

% Compute the gamma function %
g = fast_gamma(delay, dispersion, Tres, 0, TimeWindow);

% Convolve with boxcar %
if(bcw > 0)
  nbcw = round(bcw/Tres);
  b = ones(nbcw,1);
  a = conv(g,b)*Tres;
  TimeWindow = length(a)*Tres;
else
  a = g;
end

% Compute the X for an FIR %
Xfir = fast_sched2Xfir(tPres,ntrs,TR,Tres,0,TimeWindow,W);
if(isempty(Xfir)) return; end;

% Multiply the assumed shape with the FIR %
X = Xfir*a;

return;
