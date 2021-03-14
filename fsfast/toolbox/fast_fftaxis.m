function [fftaxis, deltafreq, indaxis] = fast_fftaxis(Ntp,TR)
% [fftaxis deltafreq indaxis] = fast_fftaxis(Ntp,TR)
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
%


%
% fast_fftaxis.m
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

if(nargin ~= 2) 
  msg = 'USAGE: [fftaxis deltafreq indaxis] = fast_fftaxis(Ntp,TR)';
  qoe(msg); error(msg);
end

nn = 0:round(Ntp/2);
freqmax = (1/TR)/2;         % Nyquist
deltafreq = freqmax/(Ntp/2); % Measured from 0 to Nyquist
fftaxis = deltafreq*nn;
indaxis = nn+1;
return;
