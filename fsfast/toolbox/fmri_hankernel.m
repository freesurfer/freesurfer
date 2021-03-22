function HK = fmri_hankernel(HanRadius)
%
% HK = fmri_hankernel(HanRadius)
%
% Generates a radially symetric hanning filter.
% The filter HK is a nHK x nHK matrix, where
% nHK = 2*(floor(HanRadius+.999)) + 1
% The value of the m,n th element is
%      HK(m,n) = 0.5 + cos(r*pi/2.0);
% where r is the distance of m,n from the
% center of kernel divided by the Hanning Radius.
% If r is > 1, then the component is set to zero.
% The coefficients are scaled so that the 2nd
% norm = 1.
%
%
%


%
% fmri_hankernel.m
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

if(nargin ~= 1)
  msg = 'USAGE: HK = fmri_hankernel(HanRadius)';
  qoe(msg); error(msg);
end
if(HanRadius <= 0)
  msg = 'HanRadius must be > 0';
  qoe(msg); error(msg);
end

irad = floor(HanRadius+.999);
nHK = 2*irad + 1;

HK = zeros(nHK,nHK);

for i = -irad : irad,
  m = i + irad + 1;
  for j = -irad : irad,
    n = j + irad + 1;
    r = sqrt(i*i + j*j)/HanRadius;
    if( r <= 1)
      HK(m,n) = 0.5 + cos(r*pi/2.0);
    else
      HK(m,n) = 0;
    end
  end
end

%% Normalize so that L2 = 1%%
ss = sum(reshape1d(HK).^2);
HK = HK/ss;

return;
