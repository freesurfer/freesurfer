function [fftaxis, deltafreq] = fast_fftaxis(Ntp,TR)
% [fftaxis, deltafreq] = fast_fftaxis(Ntp,TR)
%
% Returns the frequencies at which the fft is computed, 
% from DC to the nyquist frequency. There will be Ntp/2 + 1
% frequencies.
%
% Example:
% fftaxis_pos = fast_fftaxis(Nf,TR);
%
% To get the negative frequencies:
% fftaxis_neg = -fliplr(fftaxis_pos(2:end-1));
%
% $Id: fast_fftaxis.m,v 1.3 2004/01/17 05:33:39 greve Exp $

if(nargin ~= 2) 
  msg = 'USAGE: [fftaxis, deltafreq] = fast_fftaxis(Ntp,TR)';
  qoe(msg); error(msg);
end

nn = 0:round(Ntp/2);
freqmax = (1/TR)/2;         % Nyquist
deltafreq = freqmax/(Ntp/2); % Measured from 0 to Nyquist
fftaxis = deltafreq*nn;

return;
