function [SCM, nNNCond] = fmri_par2scm(Par, nTP, TR, nHEst, tPreStim, TPExclude)
%
% [SCM nNNCond] = fmri_par2scm(Par, nTP, TR, nHEst, tPreStim, <TPExclude>)
%
% Computes the Stimuli Convolution Matrix for multiple runs
%
% Arguments:
%  1. Par - this is a list of time/Simulus Id pairs (ie,
%     the Par file).  (nTPx2xnRuns)
%  2. nTP - number of time points (ie, number of samples)
%  3. TR - Time between samples 
%  4. nHEst - the number of elements of the hemodynmic impulse
%     response that will be estimated (excluding prestim).
%  5. tPreStim - time before stimulus onset to fit; this time comes
%     out of the total time window.
%  6. TPExclude - (optional) nTPxnRuns matrix with ones at locations
%     of TimePoints (ie, samples or scans) to exclude.  If TPExclude
%     does not appear as an argument, no points are excluded.
%
% Returns:
%  1. SCM - deconvolution matrix (nTP x nHEst x nRuns)
%  2. nNNCond - number of non-null (ie, non-fix) conditions.

% Notes:
%   1. If presentations do not occur simultanteously with sampling
%      (ie, the time in the parfile is not an integer multiple of
%      the TR), the presentation time is rounded to the nearest
%      sampling time.
%   2. It is assumed that the first scan event begins at t=0.
%
%


%
% fmri_par2scm0.m
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


if(nargin ~= 5 & nargin ~= 6)
  msg = 'Incorrect number of arguments.';
  qoe(msg);  error(msg);
end

nRuns = size(Par,3);

if(nargin == 5) TPExclude = zeros(nTP,2,nRuns);
else
  if(size(TPExclude,2) ~= nRuns)
    msg = 'Par and TPExclude dimensions are inconsistent';
    qoe(msg); error(msg);
  end
  if(size(TPExclude,1) ~= nTP)
    msg = 'TPExclude rows does not equal nTP';
    qoe(msg); error(msg);
  end
end

% Make sure that the number of events is not greater than the 
% number of samples.
if(size(Par,1) > nTP)
  msg = 'Number of events is greater than the number of samples';
  qoe(msg); error(msg);
end
%% Check for events before zero %%
if(length(find(Par(:,1,:)<0)) > 0)
  msg = 'Events found before t=0';
  qoe(msg); error(msg);
end
%% Check for events beyond the last sample %%
tMax = (nTP-1)*TR;
if(max(max(Par(:,1,:))) > tMax)
  msg = 'Events found beyond the last sample';
  qoe(msg); error(msg);
end

%% Compute the number of StimTypes %%
nStimTypes = max(reshape1d(Par(:,2,:))); % including fixation %
if(nStimTypes < 1)
  msg = 'There must be at least 2 Stimulus Types';
  qoe(msg); error(msg);
end

nNNCond = max(reshape1d(Par(:,2,:))) - min(reshape1d(Par(:,2,:)));

nPreStim = floor(tPreStim/TR);

for r = 1:nRuns,

  % convert the stimulus times to sample indicies
  iStimSeq = round(Par(:,1,r)/TR) + 1;

  X =  [];
  c1 = zeros(1,nHEst);
  FixationId = 0;

  for StimId = 0 : nStimTypes,

    if StimId ~= FixationId % ignore fixation

      % get the indicies of all the excitations for StimId %
      StimIndex = find( Par(:,2,r) == StimId );

      % convert these to sample indicies %
      iStim = iStimSeq(StimIndex);

      % construct the conv mtx for this Stimulus Type %
      Pulses = zeros(nTP+nPreStim,1);
      Pulses(iStim) = 1;
      c1(1) = Pulses(1);
      E = toeplitz(Pulses,c1);

      % add to global conv mtx %
      X = cat(2,X,E);

    end %%%% if StimId ~= 0 %%%

  end %%%%%% for bin  %%%%%%%

  if(nPreStim ~= 0)
    X = X(1+nPreStim:nTP+nPreStim,:);
  end

  %% Exclude specified data points %%
  iTPExclude = find(TPExclude(:,r)==1);
  if(~isempty(iTPExclude))
    fprintf(1,'   Run %d: Excluding %d Data Points:\n',r,length(iTPExclude));
    fprintf(1,'   ');
    fprintf(1,'%d ',iTPExclude);
    fprintf(1,'\n',iTPExclude);
    X(iTPExclude,:) = zeros(length(iTPExclude),size(X,2));
  end

  SCM(:,:,r) = X;

end

return;
