function indfreq = fast_indfreq(freq,Ntp,TR)
% indfreq = fast_indfreqfft(freq,Ntp,TR)
%
% Returns the index of the fft whose frequency
% is closest to freq.
%
% $Id: fast_indfreq.m,v 1.1 2003/05/22 19:39:41 greve Exp $

if(nargin ~= 3) 
  msg = 'USAGE: indfreq = indfreqfft(freq,Ntp,TR)';
  qoe(msg); error(msg);
end

freqaxis = fast_fftaxis(Ntp,TR);
[mm indfreq] = min(abs(freq-freqaxis);

return;
