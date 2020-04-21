function [F hc] = fast_lpfmtx(cutoffHz,TR,ntp,nc)
% [F hc] = fast_lpfmtx(cutoffHz,TR,ntp,<nc>)
%
% Creates a low-pass filter matrix with the given cutoff
% frequency. The filter has linear phase.
%
% hc are the filter coefficients (and the impulse response).
%
% nc is the number of coefficients, defaults to ntp. The higher nc
% is the flatter the pass-band will be and there will be more
% rejection in the stop-band, but the worse the edge effects will
% be. 
%

% Ref:
% http://www.exstrom.com/journal/sigproc

if(nargin < 3 | nargin > 4)
  fprintf('[F hc] = fast_lpfmtx(cutoffHz,TR,ntp,<nc>)\n');
  return;
end

if(~exist('nc','var')) nc = []; end
if(isempty(nc)) nc = ntp; end

wcutoff = TR*cutoffHz*2*pi;

hc = zeros(nc,1);
d1 = (nc-1)/2;
for k = 0:nc-1
  d2 = k - d1;
  if(d2 == 0) h = wcutoff / pi;
  else        h = sin( wcutoff * d2 ) / ( pi * d2 );
  end
  hc(k+1) = h;
  %fprintf('%1.15f\n', h );
end

irf = zeros(ntp,1);
irf(1:nc) = hc;
a = zeros(1,ntp);
a(1) = irf(1);
F = toeplitz(irf,a);

return;

% Testing --------------------------------------
TR = 2;
cutoffHz = .1;
ntp = 100;
nc = round(ntp/10);

[F hc] = fast_lpfmtx(cutoffHz,TR,ntp,nc);

figure(1);
plot(TR*[0:nc-1],hc);
title('Temporal Impulse Response');
xlabel('Time (sec)');

figure(2);
[fftaxis deltafreq indaxis] = fast_fftaxis(ntp,TR);
% FFT of impulse response
hcntp = zeros(ntp,1);
hcntp(1:nc) = hc;
Fhc = fft(hcntp);
% Emperical test
r = F*randn(ntp,10000);
rfft = mean(abs(fft(r)),2);
rfft = rfft(indaxis);
rfft = max(abs(Fhc(indaxis)))*rfft/max(rfft);

subplot(2,1,1);
plot(fftaxis,abs(Fhc(indaxis)),fftaxis,rfft);
title('Spectral Impulse Response');
xlabel('Frequency (Hz)');
legend('Theoretical','Emperical');
subplot(2,1,2);
plot(fftaxis,unwrap(angle(Fhc(indaxis))));






