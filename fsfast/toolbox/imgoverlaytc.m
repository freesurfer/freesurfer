function [tcimg,cmap,cscale] = imgoverlaytc(img,overlay,ovmin,ovmax,tail,rescale)
%
% [tcimg,cmap,cscale] = imgoverlaytc(img,overlay,ovmin,ovmax,tail,<rescale>)
%
% Creates a true-color compositive image with img as a gray-scale
% base and overlay as a colored overlay.
% 
% $Id: imgoverlaytc.m,v 1.1 2003/03/04 20:47:41 greve Exp $

if(nargin ~= 5 & nargin ~= 6)
  msg = 'USAGE: [imgov cmap cscale] = imgoverlaytc(img,overlay,ovmin,ovmax,tail,<rescale>)'
  qoe(msg);error(msg);
end

if(nargin == 5) rescale = 0; end

nImgRows  = size(img,1);
nImgCols  = size(img,2);

ncmap = 64;
cmgray = gray(ncmap);
cmap   = jet(ncmap);

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
ind = find(p > ovmax);
p(ind) = ovmax;

ind = find(p >= ovmin   &   p <= ovmax );

pb = p(ind);
pc = round((ncmap-1)*(pb-ovmin)/(ovmax-ovmin))+1;
tcimg(ind,:) = cmap(pc,:);

tcimg = reshape(tcimg,[size(img) 3]);

dov = (ovmax-ovmin)/(ncmap-1);
cscale = [ovmin:dov:ovmax];

return;
