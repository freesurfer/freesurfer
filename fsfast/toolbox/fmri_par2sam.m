function SAM = fmri_par2sam(Par, nTP, TR, nHEst, Baseline)
%
% SAM = Par2SAM(Par, nTP, TR, nHEst, <Baseline>)
%
% Computes the Selective Averaging Matrix
%
% Arguments:
%  1. Par - this is a list of time/SimulusId pairs (ie,
%      the Par file).
%  2. nTP - number of time points (ie, number of samples)
%  3. TR - Time between samples 
%  4. nHEst - the number of elements of the hemodynmic impulse
%     response that will be estimated.
%  5. Baseline - compensate for mean. Default is to remove
%       baseline. Set Baseline = 0 to prevent.
%
% Returns:
%  1. SAM - selective averaging matrix
%
% Notes:
%   1. If presentations do not occur simultanteously with sampling
%      (ie, the time in the parfile is not an integer multiple of
%      the TR), the presentation time is rounded to the nearest
%      sampling time.
%   2. It is assumed that the first sample occurs at t=0.
%
% Douglas N. Greve
% 1999-01-07


%
% fmri_par2sam.m
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

% Make sure that the number of events is not greater than the 
% number of samples.
if(size(Par,1) > nTP)
  error('Number of events is greater than the number of samples');
  return;
end
%% Check for events before zero %%
if(length(find(Par(1,:)<0)) > 0)
  error('Events found before t=0');
  return;
end
%% Check for events beyond the last sample %%
tMax = (nTP-1)*TR;
if(length(find(Par(:,1)>tMax)) > 0)
  error('Events found beyond the last sample');
  return;
end

%% Shift Condition Numbers to begin at 0 %%
% Par(:,2)  = Par(:,2) - min(Par(:,2));

%% Get number of stimulus types %%
nStimTypes = max(Par(:,2)); % excluding fixation %
if(nStimTypes < 1)
  error('There must be at least 1 Stimulus Types');
  return;
end

% convert the stimulus times to sample indicies
iStimSeq = round(Par(:,1)/TR) + 1;

X =  [];
c1 = zeros(1,nHEst);

for StimId = 0 : nStimTypes,

    % get the indicies of all the excitations for StimId %
    StimIndex = find( Par(:,2) == StimId );

    % convert these to sample indicies %
    iStim = iStimSeq(StimIndex);

    % construct the conv mtx for this Stimulus Type %
    Pulses = zeros(nTP,1);
    Pulses(iStim) = 1;
    c1(1) = Pulses(1);
    E = toeplitz(Pulses,c1);

    % add to global conv mtx %
    X = cat(2,X,E);

end %%%%%% for bin  %%%%%%%

%% Set Default for Baseline removal %%
if(nargin ~= 5) Baseline = 1; end

if(Baseline ~= 0)  SAM = [X ones(size(X,1),1)];
else               SAM = X;
end

return;
