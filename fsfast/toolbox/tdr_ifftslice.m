function img = tdr_ifftslice(kimg)
% img = tdr_fftslice(kimg)
%
% Performs a 2D inverse FFT of a stack of slices. The stack
% can have any number of dimensions. Does not take abs().
%
%


%
% tdr_ifftslice.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:35 $
%    $Revision: 1.3 $
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

img = [];
if(nargin ~= 1)
  fprintf('img = tdr_ifftslice(kimg)\n');
  return;
end

kimg = fftshift(fftshift(kimg,1),2);
img = ifft2(kimg);
img = fftshift(fftshift(img,1),2);


return;

