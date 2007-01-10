function indfreq = indfreqfft(freq,Ntp,TR)
% indfreq = indfreqfft(freq,Ntp,TR)


%
% fmri_indfreqfft.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
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
