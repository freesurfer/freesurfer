function img1d = imgto1d(img,keeplast)
%
% img1d = imgto1d(img,keeplast)
%
% Reshapes the image to 1D matrix. If keeplast=1,
% then the image is reshaped to a 2D image with
% the first N-1 dimension in the rows and the last
% dimension as the columns.
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
