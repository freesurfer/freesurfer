function [indpeak twfPer indtrough] = peakfinder(twf,fmin,fmax)
% [indpeak twfPer indtrough] = peakfinder(twf,<fmin>,<fmax>)
%
% Finds local peaks in a near-periodic waveform. This fails if the
% peaks are closer than 1/2 period. twfPer is the period of twf in
% samples. If the 1st or last peak are less than 90% of the mean of
% the rest of the peaks, they are exluded. 
%
% fmin,fmax are constraints on the fundamental frequency and are
% given in units of items per time point (NOT IN Hz!).
%
% $Id: peakfinder.m,v 1.9 2010/05/18 19:31:44 greve Exp $

indpeak = [];

if(nargin < 1 | nargin > 3)
  fprintf('[indpeak twfPer indtrough] = peakfinder(twf,<fmin>,<fmax>)\n');
  return;
end

Ntp = length(twf);
nn = 1:Ntp;

% detrend - necessary?
X = fast_polytrendmtx(1,Ntp,1,3);
%twf = twf - X*(inv(X'*X)*(X'*twf));

% Get major period of waveform
[fftaxis, deltafreq] = fast_fftaxis(Ntp,1);
nfft = length(fftaxis);
nnfft = 1:nfft;
twffft = abs(fft(twf-mean(twf)));
if(exist('fmin','var'))
  indok = find(fftaxis >= fmin & fftaxis <= fmax);
else
  indok = [1:length(fftaxis)];
end
[tmp k] = max(twffft(indok));
twfFreq = fftaxis(indok(k));
twfPer = 1/twfFreq;
twfPerSamp      = round(twfPer);
twfHalfPerSamp  = round(twfPerSamp/2);
twfQuartPerSamp = round(twfPerSamp/4);

%plot(fftaxis,twffft(nnfft))
%keyboard

% Assume global peak is a local peak
[tmp k0] = max(twf);
%fprintf('global peak at %d %f\n',k0,k0/25);

% Look ahead, starting at global peak
indpeak = k0;
kprev = k0;
while(1)
  % Find next max by searhing over time starting at 1/2
  % period beyond the previous max and ending one period
  % later. This fails if the peaks are closer than 1/2
  % period.
  kstart = kprev + twfHalfPerSamp;
  k = kstart + [0:twfPerSamp-1];
  indok = find(k < Ntp);
  if(length(indok) < twfQuartPerSamp) break; end
  k = k(indok);
  [tmp mmax] = max(twf(k));
  kmax = kstart + mmax - 1;
  indpeak = [indpeak kmax];
  kprev = kmax;
end

% Look behind (reverse and look ahead)
twfrev = flipud(twf(:));
k0rev = Ntp - k0 + 1; % DONT = max(twfrev);
indpeakrev = []; % dont include k0rev here
kprev = k0rev;
while(1)
  kstart = kprev + twfHalfPerSamp;
  k = kstart + [0:twfPerSamp-1];
  indok = find(k < Ntp);
  if(length(indok) < twfQuartPerSamp) break; end
  k = k(indok);
  [tmp mmax] = max(twfrev(k));
  kmax = kstart + mmax - 1;
  indpeakrev = [indpeakrev kmax];
  kprev = kmax;
end

% Convert reversed indices to forard indices
indpeakrevfor = Ntp - indpeakrev + 1;
indpeak = sort([indpeak indpeakrevfor]);
npeaks = length(indpeak);

% Decide whether to eliminate the first peak
% Compute mean of closest 3 peaks
peakfirst = twf(indpeak(1));
indpm = [2:min(4,npeaks)];
peakmean = mean(twf(indpeak(indpm)));
% Must be greater than 0.7 times this mean
if(peakfirst < .7*peakmean)  indpeak = indpeak(2:end); end

% Decide whether to eliminate the last peak
% Compute mean of closest 3 peaks
peaklast  = twf(indpeak(end));
indpm = [max(npeaks-3,1),max(npeaks-1,1)];
peakmean = mean(twf(indpeak(indpm)));
% Must be greater than 0.7 times this mean
if(peaklast  < .7*peakmean)  indpeak = indpeak(1:end-1); end

% Make sure they are unique (why not always?)
indpeak = unique(indpeak);

npeaks = length(indpeak);
indtrough = zeros(size(indpeak));
for nthpeak = 1:npeaks
  i1 = indpeak(nthpeak);
  if(nthpeak < npeaks)
    i2 = indpeak(nthpeak+1);
  else
    i2 = Ntp;
  end
  if(i2 > i1 + 1) i1 = i1 + 1; end
  [mmin imin] = min(twf(i1:i2));
  indtrough(nthpeak) = imin+i1-1;
end

if(0)
nn = 1:Ntp;
plot(nn,twf,nn(indpeak),twf(indpeak),'*',nn(indtrough),twf(indtrough),'o');
fprintf('Period %f\n',twfPerSamp);
keyboard
end

return;
