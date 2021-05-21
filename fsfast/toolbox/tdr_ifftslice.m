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

img = [];
if(nargin ~= 1)
  fprintf('img = tdr_ifftslice(kimg)\n');
  return;
end

kimg = fftshift(fftshift(kimg,1),2);
img = ifft2(kimg);
img = fftshift(fftshift(img,1),2);


return;

