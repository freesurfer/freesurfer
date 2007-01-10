function [ibrain, iair, thresh] = fmri_segment(img)
%
% [ibrain iair thresh] = fmri_segment(img)
%
% Segments the brain from the air using a single image.
% The algorithm simply measures the mean of the image
% and assigns everything over the mean to be brain
% and everything below to be air.
%
%


%
% fmri_segment.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
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
  msg = 'Usage: [ibrain iair thresh] = fmri_segment(img)';
  qoe(msg);error(msg);
end

[nRows nCols nRuns] = size(img);
nVoxels = nRows*nCols;
if(nRuns > 1) img = mean(img,3); end

thresh = mean(reshape(img,[nVoxels 1]));

ibrain = find(img >  thresh);
iair   = find(img <= thresh);

thresh = repmat(thresh, [1 nRuns]);
ibrain = repmat(ibrain, [1 nRuns]);
iair   = repmat(iair,   [1 nRuns]);

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute histogram of the image
[h v] = hist(img, nVoxels/50);

% Find the indicies where the derivative of
% the histogram is positive
iPos_dH = find(sign(diff(h))==1);

%% Compute the segmentation threshold as the first point where the
%% derivative of the histogram is positive 
if isempty(iPos_dH) thresh = max(Img1d);
else                thresh = v(iPos_dH(1));
end;

% Binarize so that subthreshold voxels are 0 and suprathreshold
% voxels are 1.
ibrain = find(img >  thresh);
iair   = find(img <= thresh);


return;
