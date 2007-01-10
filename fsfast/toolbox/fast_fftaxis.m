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
%


%
% fast_fftaxis.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
%    $Revision: 1.4 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

if(nargin ~= 2) 
  msg = 'USAGE: [fftaxis, deltafreq] = fast_fftaxis(Ntp,TR)';
  qoe(msg); error(msg);
end

nn = 0:round(Ntp/2);
freqmax = (1/TR)/2;         % Nyquist
deltafreq = freqmax/(Ntp/2); % Measured from 0 to Nyquist
fftaxis = deltafreq*nn;

return;
