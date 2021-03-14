function img1d = imgto1d(img,keeplast)
%
% img1d = imgto1d(img,keeplast)
%
% Reshapes the image to 1D matrix. If keeplast=1,
% then the image is reshaped to a 2D image with
% the first N-1 dimension in the rows and the last
% dimension as the columns.
%


%
% imgto1d.m
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

img1d = [];

if(nargin ~= 1 & nargin ~= 2)
  fprintf('USAGE: img1d = imgto1d(img,<keeplast>)\n');
  return;
end

if(nargin == 1) keeplast = 0; end

sz = size(img);
nsz = length(sz);

if(~ keeplast)
  newsz = [prod(sz) 1];
else
  newsz = [prod(sz(1:nsz-1)) sz(nsz)];
end

img1d = reshape(img, newsz);


return;
