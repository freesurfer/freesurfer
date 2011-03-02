function indfreq = indfreqfft(freq,Ntp,TR)
% indfreq = indfreqfft(freq,Ntp,TR)


%
% fmri_indfreqfft.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:06 $
%    $Revision: 1.5 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

if(nargin ~= 3) 
  msg = 'USAGE: indfreq = indfreqfft(freq,Ntp,TR)';
  qoe(msg); error(msg);
end

df = 1/(Ntp*TR); % frequency increment

fr = freq/df + 1; % ratio of frequency to increment

if( rem(fr,1) ~= 0 )
  % indfreq = [floor(fr) ceil(fr)];
  indfreq = round(fr);
else
  indfreq = fr;
end

return;
