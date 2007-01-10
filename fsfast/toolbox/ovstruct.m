function [ovs] = ovstruct(img,actimg,ovmin,ovmax)
%
% [ovs] = ovstruct(img,actimg,ovmin,ovmax)
%
%


%
% ovstruct.m
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

if(nargin ~= 4)
  msg = 'USAGE: [ovs] = ovstruct(img,actimg,ovmin,ovmax)';
  qoe(msg);error(msg);
end

nImgRows  = size(img,1);
nImgCols  = size(img,2);

ovs.min = ovmin;
ovs.max = ovmax;
ovs.nframes = size(actimg,3);

ncmap = 64;
cmap  = jet(ncmap);

actimg = reshape1d(fmri_imresize(actimg, [nImgRows nImgCols]));

for frame = 1:nframes,
  p = actimg(:,:,frame);

  % Saturate at ovmax %
  ovs.pos[frame].ind = find(p > ovmax);
  p(ind) = ovmax;

ind = find(p >= ovmin   &   p <= ovmax );

pb = p(ind);
pc = round((ncmap-1)*(pb-ovmin)/(ovmax-ovmin))+1;
tcimg(ind,:) = cmap(pc,:);

tcimg = reshape(tcimg,[size(img) 3]);

dov = (ovmax-ovmin)/(ncmap-1);
cscale = [ovmin:dov:ovmax];

return;
