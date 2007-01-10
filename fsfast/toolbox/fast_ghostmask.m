function ghostmask = fast_ghostmask(mask)
% ghostmask = fast_ghostmask(mask)
% 
% Create a mask of the ghost from the mask of the main image
%   Each slice is an epi slice
%   The phase encode changes from one row to another (beware
%     matlabs row major!)
% mask is of size nrows X ncols X nslices (as is ghostmask)
% mask does not have to be a binary mask. The voxel values can 
%   be anything.
% 
% To compute the row that a given brain voxel maps to in the ghost:
%   if(brainrow <= N/2) ghostrow = brainrow + N/2
%   else                ghostrow = brainrow - N/2
%   end
%
%


%
% fast_ghostmask.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
%    $Revision: 1.3 $
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

ghostmask = [];

if(nargin ~= 1)
  fprintf('ghostmask = fast_ghostmask(mask)\n');
  return;
end

[nrows ncols nslices] = size(mask);

nrows2 = round(nrows/2); % half the number of rows
indA = 1:nrows2;         % indices of the top half
indB = nrows2+1:nrows;   % indices of the bottom half
    
ghostmask = mask([indB indA],:,:);

return;
