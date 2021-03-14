function h = imagesc1d(img1d,nRows,nCols)
%
% h = imagesc1d(img1d,nRows,nCols)
%
% Produces a 2-D image of a 1-D vector.  If nRows or nCols
% are unspecified, it attempts to generate a square image.
%
%


%
% imagesc1d.m
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
