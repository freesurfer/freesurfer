function h = imagesc1d(img1d,nRows,nCols)
%
% h = imagesc1d(img1d,nRows,nCols)
%
% Produces a 2-D image of a 1-D vector.  If nRows or nCols
% are unspecified, it attempts to generate a square image.
%
% $Id: imagesc1d.m,v 1.1 2003/03/04 20:47:41 greve Exp $


% number of pixels in the 1-D image
nPixels = length(img1d);

% check whether rows and cols were specified
if(nargin < 3)
  nSide = floor(sqrt(nPixels));
  nRows = nSide;
  nCols = nSide;
  if(nSide*nSide ~= nPixels)
    error('Length of Image must be square or spec rows and cols');
    return;
  end
  else
  if(nRows*nCols ~= nPixels)
    error('Length of Image must be = nRows * nCols');
    return;
  end
end


h = imagesc(reshape(img1d, [nRows nCols]));

return;
