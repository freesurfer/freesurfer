function tseq = fast_schedseq(seq, tpercond, TER, TScan, TPreScan)
%
% tseq = fast_schedseq(seq, tpercond, TER, TScan, TPreScan)
%
% Computes a (random) time at which each event in sequence seq 
% can occur. No optimization is performed here.
%
% seq is the sequence of events. Each event is coded 1-N. Zero cannot
% be used as an event code. Event codes should be contiguous from 1-N.
%
% tpercond is the amount of time allocated for each event type. Its
% length must be equal to the number of event types.
%
% TER is the temporal estimation resolution (in seconds).
%
% TScan is the length of the scan (in seconds) during which data
% are collected AND stored (do not include discarded acquisitions).
%
% TPreScan is the amount of time before the onset of scanning that
% stimuli should begin to appear. Where "the onset of scanning" is
% defined as the first image collected AND stored. Note that TPreScan
% is not inherently related to the number discarded acquisitions. 
% TPreScan allows the user to present stimuli before images are
% collected. If this is done, then these presentations should be
% considered during the analysis.
%
% tseq is a list of times at which the items in seq should be presented.
% The times are in seconds relative to the first image that is
% collected AND stored.  
%
%


%
% fast_schedseq.m
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

tseq = [];

if(nargin ~= 5)
  msg = 'USAGE: tseq = fast_schedseq(seq, tpercond, TER, TScan, TPreScan)';
  fprintf('%s\n',msg);
  return;
end

nevents = length(seq);

% Check for 0s in seq %
n = length(find(seq==0));
if(n ~= 0) 
  fprintf('ERROR (fast_schedseq): sequence contains zeros\n');
  return;
end

% Get the number of event types %
neventtypes = length(unique(seq));
if(neventtypes ~= max(seq))
  fprintf('ERROR (fast_schedseq): sequence has missing ids\n');
  return;
end

if(length(tpercond) ~= neventtypes)
  fprintf('ERROR (fast_schedseq): dimension of tpercond (%d)\n',...
          length(tpercond));
  fprintf('does not match the number of event types (%d)\n',neventtypes);
  return;
end

% Get the number of presentations of each type %
for n = 1:neventtypes
  npercond(n) = length(find(seq==n));
end

% Compute the total amount of scan time
TTot = TPreScan + TScan;

% Compute the total number ters
nters = floor(TTot/TER);

% Compute the total amount of stimulation time
TTotStim = sum(npercond .* tpercond);

% Compute the total amount of null time %
TTotNull = TTot - TTotStim;

if(TTotNull <= 0)
  fprintf('ERROR (fast_schedseq): Not enough null time. Decrease\n');
  fprintf('the number of events or the duration of the events or\n');
  fprintf('increase the amount of scan time.\n');
  return;
end

nnulls = round(TTotNull/TER);

% A list of the duration of each event in the sequnce %
dseq = tpercond(seq);

% The null time to put before each event
tpernull = TER*fast_npernull(nevents,nnulls);

% The duration of each event, incl the following null
dtseq = tpernull(2:nevents+1) + dseq(1:nevents) ;

tseq(1) = tpernull(1);
tseq(2:nevents) = tpernull(1) + cumsum(dtseq(1:nevents-1));
tseq = tseq - TPreScan;

return;
