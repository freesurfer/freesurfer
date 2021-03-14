function [tcimg,cmap,cscale] = imgoverlaytc(img,overlay,ovmin,ovmax,tail,rescale)
%
% [tcimg,cmap,cscale] = imgoverlaytc(img,overlay,ovmin,ovmax,tail,<rescale>)
%
% Creates a true-color compositive image with img as a gray-scale
% base and overlay as a colored overlay.
% 
%


%
% imgoverlaytc.m
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

if(nargin ~= 5 & nargin ~= 6)
  msg = 'USAGE: [imgov cmap cscale] = imgoverlaytc(img,overlay,ovmin,ovmax,tail,<rescale>)'
  qoe(msg);error(msg);
end

if(nargin == 5) rescale = 0; end

nImgRows  = size(img,1);
nImgCols  = size(img,2);

ncmap = 64;
cmgray = gray(ncmap);
%cmap   = jet(ncmap);
%cmap   = hot(ncmap);
%cmap    = nmrcolormap(ncmap,tail);
cmap    = nmrcolormap(ncmap,'both');

if(rescale)
  t1min = min(reshape1d(img));
  t1max = max(reshape1d(img));
  t1b   = round((ncmap-1)*(img-t1min)/(t1max-t1min))+1;
  t1c   = reshape1d(t1b);
  tcimg = cmgray(t1c,:);
else
  tcimg = cmgray(reshape1d(img),:);
end

p = reshape1d(fmri_imresize(overlay, [nImgRows nImgCols]));

switch(tail)
  case 'abs',  p = abs(p);
  case 'neg',  p = -p;
end

% Saturate at ovmax %
ind = find(abs(p) > ovmax);
p(ind) = ovmax.*sign(p(ind));

ind = find(p >= ovmin   &   p <= ovmax );

pb = p(ind);
pc = round((ncmap-1)*(pb-ovmin)/(ovmax-ovmin))+1;
tcimg(ind,:) = cmap(pc,:);

tcimg = reshape(tcimg,[size(img) 3]);

dov = (ovmax-ovmin)/(ncmap-1);
cscale = [ovmin:dov:ovmax];

return;
