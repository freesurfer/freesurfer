function [bpar,nstimlen] = parboxcar(erpar,TER,stimlen,runlen)
% [bpar nstimlen] = parboxcar(erpar,TER,<stimlen>,<runlen>)
%
% Convolves stimulus onset times with boxcar in order to represent
% non-infinitesimal stimulus lengths. 
%
% A single stimulus presentation is replaced with a series of
% presentations separated in time by TER and spanning the duration of
% the stimulus. The stimulus duration can be specified in one of three
% ways. (1) all stimuli have the same length specified by stimlen. (2)
% each stimulus can have its own length specified by the 3rd column of
% erpar. (3) each stimulus can have its own length specified by the
% time between its onset and the onset of the next stimulus. For the
% 3rd option, the duration of the last stimulus is ambiguous if that
% stimulus is non-null. In this case, add a null stimulus at the time
% of the end of the last stimulus or specify the length of the run
% with run len.
%
% erpar - stimulus timing, nstim-by-2 or nstim-by-3:  
%   1: times of stimulus onsets
%   2: stimulus ids
%   3: stimulus duration (or use stimlen)
% TER - in seconds
% stimlen - length of stimulus in seconds if not done
%   per stimulus. When per-stim duration and runlen need
%   to be specified, use simtlen=[].
% runlen - duration of the run in sec. May be needed when
%   automatically determining the duration of each stimulus from
%   the stimulus onset times.
% 
% bpar - nbstim-by-2.
% nstimlen - number TERs for each stimulus.
%
% Notes:
%  1. Stim Id 0 in erpar will not be represented in the output. 
%
%
%


%
% parboxcar.m
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

bpar = [];
nstimlen = [];
if(nargin < 1 | nargin > 4)
  fprintf('bpar = parboxcar(erpar,TER,<stimlen>,<runlen>)\n');
  return;
end
if(~exist('stimlen','var')) stimlen=[]; end
if(~exist('runlen','var'))  runlen=[]; end

erpar0 = erpar;

nstim = size(erpar,1);

% Get a list of stimulus durations, either per-stimulus
% from erpar(:,3), or apply the scalar stimlen to all stims 
if(~isempty(stimlen))
  % Use same stim len for all stimuli
  nstimlen = round(stimlen/TER)*ones(nstim,1);
elseif(size(erpar,2) == 3)
  % Use per-stimulus length
  nstimlen = round(erpar(:,3)/TER);
else
  % Auto extract per-stimulus lengths from stim timing
  if(~isempty(runlen))
    % Remove any stimuli presented beyond the end of the run
    indkeep = find(erpar(:,1)<runlen);
    erpar = erpar(indkeep,:);
    nstim = size(erpar,1);
    % Add null stimulus at the end to signal end-of-run
    erpar(nstim+1,1) = runlen;
    erpar(nstim+1,2) = 0;
    nstim = size(erpar,1);
  end
  if(erpar(end,2) ~= 0)
    fprintf('ERROR: last condition must be 0 for auto.\n');
    return;
  end
  nstimlen = zeros(nstim-1,1);
  for nthstim = 1:nstim-1
    dt = erpar(nthstim+1,1)-erpar(nthstim,1);
    nstimlen(nthstim) = round(dt/TER);
    if(nthstim ~= nstim-1 & nstimlen(nthstim)==0)
      fprintf('ERROR: parboxcar(): stimulus %d at time %g has duration %g < TER=%g\n',...
	      nthstim,erpar(nthstim,1),dt,TER);
      keyboard
      return;
    end
  end
end

nthbpar = 1;
for nthstim = 1:nstim
  id = erpar(nthstim,2);
  if(id == 0) continue; end

  t  = erpar(nthstim,1);
  for nthter = 1:nstimlen(nthstim)
    bpar(nthbpar,1) = t;
    bpar(nthbpar,2) = id;
    t = t + TER;
    nthbpar = nthbpar + 1;
  end
  
end

return;


















