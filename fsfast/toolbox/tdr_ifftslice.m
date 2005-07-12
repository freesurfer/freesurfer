function img = tdr_ifftslice(kimg)
% img = tdr_fftslice(kimg)
%
% Performs a 2D inverse FFT of a stack of slices. The stack
% can have any number of dimensions. Does not take abs().
%
% $Id: tdr_ifftslice.m,v 1.2 2005/07/12 19:02:47 greve Exp $

img = [];
if(nargin ~= 1)
  fprintf('img = tdr_ifftslice(kimg)\n');
  return;
end

kimg = fftshift(fftshift(kimg,1),2);
img = ifft2(kimg);
img = fftshift(fftshift(img,1),2);


return;

