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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
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
