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
