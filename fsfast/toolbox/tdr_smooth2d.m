function [imgsm, G] = tdr_smooth2d(img,sig1,sig2)
% [imgsm filtimg] = tdr_smooth2d(img,sig1,sig2)
% 
% 2D gaussian smoother.
%
% sig1 - gaussian standard dev for the rows
% sig2 - gaussian standard dev for the columns
%
% $Id: tdr_smooth2d.m,v 1.1 2003/11/06 19:45:37 greve Exp $

if(nargin ~= 3)
  fprintf('[imgsm filtimg] = tdr_smooth2d(img,sig1,sig2)\n');
  return;
end

[nrows ncols] = size(img);

% Create a 2d kernel - same in image and kspace 
%filtimg = exp(-1/2*(((([1:nrows]-(nrows/2+1))'*ones(1,ncols))/nrows*sig1).^2+((ones(nrows,1)*([1:ncols]-(ncols/2+1)))/ncols*sig2).^2));

x = ([1:nrows] - (nrows/2 + 1))/sig1;
xx = ones(nrows,1) * x;
y = ([1:ncols]' - (ncols/2 + 1))/sig2;
yy = y * ones(1,ncols);
G = exp( -(xx.^2)/2 + -(yy.^2)/2 );
G = G/max(reshape1d(G));

%sfiltimg = fftshift(fftn(ifftshift(filtimg)));
%figure; imagesc(log(filtimg)); colorbar; axis equal; axis image;
%figure; imagesc(log(abs(sfiltimg))); colorbar; axis equal; axis image;

% Conver the image to kspace
%kimg = fftshift(ifftn(img));

% Smooth in kspace and return to image space
%imgsm = fftn(ifftshift(kimg.*filtimg));

kimg  = fftn(img);
kG    = fftn(G);
imgsm = fftshift(ifftn(kimg.*kG));

if( max(abs(reshape1d(imag(img)))) == 0 )
  % Only take the magnitude if input is totally really
  imgsm = abs(imgsm);
end

%keyboard

