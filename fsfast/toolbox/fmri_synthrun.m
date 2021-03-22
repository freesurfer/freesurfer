function [srs, fSlices] = fmri_synthrun(srs)
%
%


%
% fmri_synthrun.m
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

srs.Seed = sum(100*clock);
randn('state',srs.Seed);
rand('state',srs.Seed);
fprintf('Seed = %d\n',srs.Seed);

[Mss srs.Rss]   = fmri_subsampmat(srs.TR,srs.Ntp,srs.TER);
srs.Nts = srs.Ntp * srs.Rss;
Nh = floor(srs.TimeWindow/srs.TER);

if(isempty(srs.Par))
  Nc = length(srs.NPerCond) + 1;
  srs.Nswap = round(srs.Nts/2);
  N0 = srs.Nts - sum(srs.NPerCond);
  seqbase = fmri_seqbase([N0 srs.NPerCond]);
  fprintf('\nSearching for Optimal Seq over %d\n',srs.Nsearch);
  [seq tr travg trstd] = fmri_minseqtr(seqbase,Nh,srs.Nsearch,1,srs.Nswap,Mss);
  fprintf('Min Trace = %g\n',tr);
  srs.Par = [ srs.TER*[0:srs.Nts-1]' seq]; %'
else
  tmp = reshape1d(srs.Par(:,2,:));
  Nc = max(tmp) + 1;
  N0 = length(find(tmp==0));
end

Nv = srs.Nrows*srs.Ncols;
Nnnc = Nc - 1;

nRnnMaxLag = round(srs.RnnMaxDelay/srs.TER);
Rnn = arCorrFun(srs.alpha,srs.rho,nRnnMaxLag,2);

%%%% Ideal HDR %%%%%%%
t = srs.TER*[0:Nh-1]';  %'
SigAmp = srs.Offset*srs.PercentSigChange/100;
hideal = SigAmp*fmri_hemodyn(t,srs.delta,srs.tau);
hideal2 = repmat(hideal, [Nnnc 1]);

% Weight the amplitude of each condition %
wcond = [1:Nnnc]; 
wcond2 = repmat(wcond, [length(hideal) 1]);
wcond3 = reshape1d(wcond2);
%hideal2 = hideal2 .* wcond3;

SignalMask = zeros(srs.Nrows,srs.Ncols);
if(~isempty(srs.ROI))
  Rroi       = [srs.ROI(1):srs.ROI(3)];
  Croi       = [srs.ROI(2):srs.ROI(4)];
else
  Rroi = [];
  Croi = [];
  if(srs.SNR ~= 0)
    fprintf('WARNING: no ROI specified, reseting SNR to 0\n');
  end
  srs.SNR = 0;
end

SignalMask(Rroi,Croi) = 1;
SignalMask = reshape(SignalMask,[1 Nv]);
iSignal    = find(SignalMask==1);
iNoSignal  = find(SignalMask==0);
%MeanMask = (.05 + .95*SignalMask)*srs.Offset;
MeanMask = SignalMask*srs.Offset + ~SignalMask*srs.Offset/30;
MeanMask = repmat(MeanMask, [srs.Ntp 1]);

Trend = srs.Trend * srs.TR * [0:srs.Ntp-1]'; %'
Trend2 = repmat(Trend, [1 Nv]);

h  = hideal2 * SignalMask;
X0 = fmri_par2scm(srs.Par, -1, srs.Nts, srs.TER, Nh, 0);
X  = Mss*X0;
srs.traceiXtX = trace(inv(X'*X)); %'

fSignal   = X * h;
if(srs.SNR ~= 0)
  sig       = fSignal(:,iSignal(1));
  fSigVar   = cov(sig);
  fSigMean  = mean(sig);
  fSigMin   = min(sig);
  fSigMax   = max(sig);
  srs.Sig   = sig;
  srs.SigMean = fSigMean;
  srs.SigVar  = fSigVar;
  srs.SigSS   = sum(sig.^2);
  fprintf('  Sig Mean = %g, Std = %g, Var = %g, Min = %g, Max = %g, Range = %g\n',...
          fSigMean,sqrt(fSigVar),fSigVar,fSigMin,fSigMax,fSigMax-fSigMin);
  srs.SigVar = fSigVar;
  avghideal = mean(hideal2);
  scfact    = Nh * sum(srs.NPerCond)/(N0 + sum(srs.NPerCond));
  fprintf('  Avg hIdeal: %g, scaled = %g\n',avghideal,avghideal*scfact);

  exps2 = (hideal2' * X' * X * hideal2)/size(X,1); %
  triXtX = trace(inv(X'*X)); %'
  expSigVar = (exps2-(avghideal*scfact).^2);
  expErrVar = expSigVar*triXtX/srs.SNR;
end

fSlices = zeros(srs.Ntp,Nv,srs.Nslices);

if(srs.SNR > 0)  fNoiseVar = fSigVar/srs.SNR;
else             fNoiseVar = 30;
end
fprintf('  Noise Var = %g, tr(inv(XtX)) = %g, E(e2) = %g\n',...
	fNoiseVar, srs.traceiXtX,fNoiseVar*srs.traceiXtX);

for slice = 1:srs.Nslices,

  fprintf('Synthesizing slice %d ----- \n',slice);

  if(srs.SNR >= 0)

    fNoiseStd = sqrt(fNoiseVar);
    if(srs.RnnMaxDelay==0)
      fNoise  = fNoiseStd*randn(srs.Ntp,Nv);
    else
      fNoise  = fNoiseStd*fmri_randc(srs.Ntp*Nv,Rnn);
      fNoise  = reshape(fNoise,[srs.Ntp Nv]);
    end

  else
    fNoise = zeros(size(fSignal));
  end

  if(srs.SNR ~= 0)  
    fSigNoise = fSignal + fNoise;
  else
    fSigNoise = fNoise;
  end
  fSigNoise = fSigNoise + MeanMask + Trend2;

  %% Clip and Integerize %%
  indltz = find(fSigNoise<0);
  fprintf(1,'    Clipping %4d voxels to 0 \n',length(indltz));
  fSigNoise(indltz) = 0;

  indltz = find(fSigNoise>32000);
  fprintf(1,'    Clipping %4d voxels to 32000 \n',length(indltz));
  fSigNoise(indltz) = 32000;
  fSigNoise = floor(fSigNoise);

  fSlices(:,:,slice) = fSigNoise;

end

fSlices = permute(fSlices, [2 1 3]);
fSlices = reshape(fSlices, [srs.Nrows srs.Ncols srs.Ntp srs.Nslices]);

% Corrupt the first Nskip %
if(srs.Nskip > 0)
  fSlices(:,:,[1:srs.Nskip],:) = 0;
end


return;

