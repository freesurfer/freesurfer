% plthankernels.m
% Plots true hanning window as well as that used by selavg
% and selxavg


%
% plthankernels.m
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
