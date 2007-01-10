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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
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
