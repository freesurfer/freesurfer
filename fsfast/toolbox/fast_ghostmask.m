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
