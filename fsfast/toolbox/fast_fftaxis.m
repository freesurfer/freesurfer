function [fftaxis, deltafreq] = fast_fftaxis(Ntp,TR)
% [fftaxis, deltafreq] = fast_fftaxis(Ntp,TR)
%
% Returns the frequencies at which the fft is computed, 
% for DC to the nyquist frequency.
%
% $Id: fast_fftaxis.m,v 1.2 2003/12/02 23:07:08 greve Exp $

if(nargin ~= 2) 
  msg = 'USAGE: [fftaxis, deltafreq] = fast_fftaxis(Ntp,TR)';
  qoe(msg); error(msg);
end

nn = 1:round(Ntp/2);
freqmax = (1/TR)/2;         % Nyquist
deltafreq = freqmax/(Ntp/2);  % Measured from 0 to Nyquist
fftaxis = deltafreq*(nn-1);

return;
