% plthankernels.m
% Plots true hanning window as well as that used by selavg
% and selxavg


%
% plthankernels.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
%    $Revision: 1.2 $
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

r = [0:.01:1];

htrue = 1 + cos(r*pi);
htrue = htrue/max(htrue);

hselavg = 0.5 + cos(r*pi/2);
hselavg = hselavg/max(hselavg);
ind = find(hselavg < 0);
hselavg(ind) = 0;

% For Selxavg, turns out to be same as selavg
nHK = 10;
HK = fmri_hankernel(nHK);
HK2 = HK/sum(reshape1d(HK));
hk2 = HK2(nHK+1,nHK+1:2*nHK+1);
hk2 = hk2/max(hk2);
rhk2 = [0:nHK]/nHK;

plot(r,htrue,r,hselavg,rhk2,hk2);
xlabel('r/hanrad');
