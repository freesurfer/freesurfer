function img2 = fast_pixrep(img)
%
% img2 = fast_pixrep(img)
%
% Double image size with pixel replication of the rows and
% columns. Works with multiple frames.
%


%
% fast_pixrep.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
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

if(nargin ~= 1)
   msg = 'USAGE: img2 = fast_pixrep(img)';
   qoe(msg);error(msg);
end

[nrows ncols nother] = size(img);
nv = nrows*ncols;

% Create index map - the value at each pixel will be the
% index from the original image
img2 = 1:nv;
img2 = repmat(img2, [2 1]);
img2 = reshape(img2,[2*nrows ncols]);
img2 = reshape(img2', [1 2*nv]); %'
img2 = repmat(img2, [2 1]); 
img2 = reshape(img2, [2*nrows 2*ncols])'; %'

% Reshape the original image to be pixels-by-frames
img = reshape(img, [nv nother]);

% Map and Reshape
img2 = reshape(img(reshape1d(img2),:), [2*nrows 2*ncols nother]);

return;
