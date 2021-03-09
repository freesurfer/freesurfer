function [imgsm, G] = tdr_smooth2d(img,sig1,sig2)
% [imgsm filtimg] = tdr_smooth2d(img,sig1,sig2)
% 
% 2D gaussian smoother.
%
% sig1 - gaussian standard dev for the rows
% sig2 - gaussian standard dev for the columns
%
%


%
% tdr_smooth2d.m
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

if(nargin ~= 3)
  fprintf('[imgsm filtimg] = tdr_smooth2d(img,sig1,sig2)\n');
  return;
end

[nrows ncols] = size(img);
nrows2 = 2*nrows;
ncols2 = 2*ncols;

x = ([1:nrows2] - (nrows2/2 + 1))/sig1;
xx = ones(nrows2,1) * x;
y = ([1:ncols2]' - (ncols2/2 + 1))/sig2;
yy = y * ones(1,ncols2);
G = exp( -(xx.^2)/2 + -(yy.^2)/2 );
G = G/sum(reshape1d(G));

% zero pad
zpimg = zeros(2*size(img));
zpimg(1:nrows,1:ncols) = img;

kimg  = fftn(zpimg); 
kG    = fftn(G);
imgsm = fftshift(ifftn(kimg.*kG));
imgsm = imgsm(1:nrows,1:ncols); 

if( max(abs(reshape1d(imag(img)))) == 0 )
  % Only take the magnitude if input is totally really
  imgsm = abs(imgsm);
end

%keyboard

