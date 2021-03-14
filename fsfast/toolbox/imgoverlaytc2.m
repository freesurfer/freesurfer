function [tcimg,cmap,cscale] = imgoverlaytc2(tcimg,overlay,ovmin,ovmax,tail,interp)
%
% [tcimg,cmap,cscale] = imgoverlaytc(imgtc,overlay,ovmin,ovmax,tail,interp)
%
% Creates a true-color composite image with img as a gray-scale
% base and overlay as a colored overlay.  Same as imgoverlaytc() but
% the input image (tcimg) is already in true-color format.
% 
%


%
% imgoverlaytc2.m
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

if(nargin ~= 6)
  msg = 'USAGE: [imgtc cmap cscale] = imgoverlaytc2(imgtc,overlay,ovmin,ovmax,tail,interp)'
  qoe(msg);error(msg);
end

ncmap = 128;
%cmap   = jet(ncmap);
%cmap   = hot(ncmap);
cmap    = nmrcolormap(ncmap,tail);

nOvRows = size(overlay,1);
nOvCols = size(overlay,2);

nImgRows  = size(tcimg,1);
nImgCols  = size(tcimg,2);
if(exist('imresize') & interp & nImgRows == nOvRows & nImgCols == nOvCols)
  % Quadrupal size of base image so can interpolate %
  tcimg = fast_pixrep(tcimg);
  tcimg = fast_pixrep(tcimg);
end

nImgRows  = size(tcimg,1);
nImgCols  = size(tcimg,2);
nv = nImgRows*nImgCols;
tcimg = reshape(tcimg, [nv 3]);

if(exist('imresize') & interp)
  p = reshape1d(imresize(overlay, [nImgRows nImgCols],'bicubic'));
else
  p = reshape1d(fmri_imresize(overlay, [nImgRows nImgCols]));
end

if(strcmp(tail,'neg') | strcmp(tail,'abs') | strcmp(tail,'pos'))

  switch(tail)
    case 'abs',  p = abs(p);
    case 'neg',  p = -p;
  end

  % Saturate at ovmax %
  ind = find(p > ovmax);
  p(ind) = ovmax;

  ind = find(p >= ovmin);

  pb = p(ind);
  pc = round((ncmap-1)*(pb-ovmin)/(ovmax-ovmin))+1;
  tcimg(ind,:) = cmap(pc,:);

  tcimg = reshape(tcimg,[nImgRows nImgCols 3]);

  dov = (ovmax-ovmin)/(ncmap-1);
  cscale = [ovmin:dov:ovmax];

elseif(strcmp(tail,'both') | strcmp(tail,'posneg'))

  % Saturate at ovmax, pos tail %
  ind = find(p > ovmax);
  p(ind) = ovmax;

  % Saturate at ovmax, neg tail %
  ind = find(-p > ovmax);
  p(ind) = -ovmax;

  % map positive tails %
  ind = find(p >= ovmin);
  pb = p(ind);
  pc = ncmap/2 + round((ncmap/2-1)*(pb-ovmin)/(ovmax-ovmin))+1;
  tcimg(ind,:) = cmap(pc,:);

  % map negative tails %
  ind = find(p <= -ovmin);
  pb = p(ind);
  pc = ncmap/2 - round((ncmap/2-1)*(-pb-ovmin)/(ovmax-ovmin))+1;
  tcimg(ind,:) = cmap(pc,:);

  tcimg = reshape(tcimg,[nImgRows nImgCols 3]);

  dov = (ovmax-ovmin)/(ncmap/2-1);
  pos_cscale = [ovmin:dov:ovmax];
  neg_cscale = [ovmax:-dov:ovmin];
  cscale = [-neg_cscale pos_cscale];

end

return;
