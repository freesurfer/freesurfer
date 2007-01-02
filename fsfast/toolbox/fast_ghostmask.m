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
% $Id: fast_ghostmask.m,v 1.2 2007/01/02 22:48:38 greve Exp $

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
