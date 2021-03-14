function [SCM, nNNCond, nPerCond] = fmri_par2scm(Par, Nconds, nTP, TR, nHEst, tPreStim, TPExclude)
%
% [SCM nNNCond nPerCond] = fmri_par2scm(Par, Nconds, nTP, TR, nHEst, tPreStim,
%                              <TPExclude>)
%
% Computes the Stimuli Convolution Matrix for multiple runs
%
% Arguments:
%  1. Par - this is a list of time/Simulus Id pairs (ie,
%     the Par file).  (nTPx2xnRuns)
%  2. Nconds - number of conditions, including null condition.
%  3. nTP - number of time points (ie, number of samples)
%  4. TR - Time between samples 
%  5. nHEst - the number of elements of the hemodynmic impulse
%     response that will be estimated (excluding prestim).
%  6. tPreStim - time before stimulus onset to fit; this time comes
%     out of the total time window.
%  7. TPExclude - (optional) nTPxnRuns matrix with ones at locations
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
%   3. Null condition is assumed to be zero.
%
%


%
% fmri_par2scm2.m
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

if(nargin ~= 6 & nargin ~= 7)
  msg = 'USAGE: fmri_par2scm(Par, Nconds, nTP, TR, nHEst, tPreStim, <TPExclude>)';
  qoe(msg);  error(msg);
end

nRuns = size(Par,3);
NullCondId = 0;

if(nargin == 6) TPExclude = zeros(nTP,2,nRuns);
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

%% Check for events before zero %%
ilz = find(Par(:,1,:)<0);
if(~isempty(ilz))
  fprintf(1,'INFO: Events found before t=0, ignoring\n');
end

%% Check for events beyond the last sample %%
tMax = (nTP-1)*TR;
if(max(max(Par(:,1,:))) > tMax)
  fprintf(1,'WARNING: Stimulus presentations found beyond the last sample... ignoring\n');
  fprintf(1,'         tParMax = %g > tSampMax =%g\n',max(max(Par(:,1,:))),tMax);
end

%% Compute the number of StimTypes %%
FirstCondId = min(reshape1d(Par(:,2,:)));
LastCondId  = max(reshape1d(Par(:,2,:))); 
if(Nconds < 0)
  Nconds = LastCondId - FirstCondId + 1;
end
nNNCond = Nconds -1;

% Get nPerCond %
nPerCond = [];
for n = 1:nNNCond + 1,
  cid = (n-1) + FirstCondId;
  nPerCond(n) = length(squeeze(find(Par(:,2,:)==cid)));
end
fprintf(1,'FirstCondId = %d, LastCondId = %d, Nconds = %d, nNNC = %d\n',...
           FirstCondId, LastCondId, Nconds ,nNNCond);
fprintf(1,'nPerCond:  ');
fprintf(1,'%5d ',nPerCond);
fprintf(1,'\n');

nPreStim = floor(tPreStim/TR);
tTR = [0:TR:tMax];


for r = 1:nRuns,

  X =  [];
  c1 = zeros(1,nHEst);

  for StimId = 1:Nconds-1,

    if StimId ~= NullCondId,

      % get the indicies of all the excitations for StimId %
      iStimIndex = find( Par(:,2,r) == StimId );

      % convert the stimulus times to sample indicies
      tPres = Par(iStimIndex,1,r);
      iok = find(tPres >= 0 & tPres <= tMax);
      tPres = tPres(iok);
      nPres = hist(tPres,tTR); 

      % construct the conv mtx for this Stimulus Type %
      Pulses = [nPres zeros(1,nPreStim) ];
      c1(1) = Pulses(1);
      E = toeplitz(Pulses,c1);

      % add to global conv mtx %
      X = cat(2,X,E);

    end %%%% if StimId ~= 0 %%%

  end %%%%%% for bin  %%%%%%%

  if(nPreStim ~= 0)
    %size(X)
    %nPreStim
    %nTP
    X = X(1+nPreStim:nTP+nPreStim,:);
    %X = X(1:nTP,:);
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
