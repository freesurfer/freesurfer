function [indm,r,c,s] = sub2indmask(r,c,s,volmaskind,indmask)
% [indm r c s] = sub2indmask(r,c,s,volmaskind,indmask)
%
% indm are the indices into the mask that correspond
% to voxels at [r c s]. 
%
% indmask = find(mask);
% volmaskind = zeros(size(mask));
% volmaskind(indmask) = 1:length(indmask);
%
% If rcs is not in the mask, then it is excluded. In this case
% the number of components in indm will be less than that in
% r,c,s. If any of the rcs are out of the volume, then they
% are also excluded. The output r, c, and s will reflect these 
% exclusions.
%
% Example: let rcs = [20 30 17] and let this voxel be in the
% mask. If there is only one other voxel in the mask at [1 1 1],
% then this function will return indm=2 since [20 30 17] will
% be the 2nd voxel in the mask.
%
%


%
% sub2indmask.m
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


indm = [];

if(nargin ~= 5)
  fprintf('[indm r c s] = sub2indmask(r,c,s,volmaskind,indmask)\n');
  return;
end

volsize = size(volmaskind);

% Remove crs that are out-of-volume
ind_invol = find( (r > 0 & r <= volsize(1)) & ...
		  (c > 0 & c <= volsize(2)) & ...
		  (s > 0 & s <= volsize(3)));

if(isempty(ind_invol)) return; end

r = r(ind_invol);
c = c(ind_invol);
s = s(ind_invol);

indvol = sub2ind(volsize,r,c,s);
indm = volmaskind(indvol);
tmp = find(indm ~= 0);
indm = indm(tmp);



return;





