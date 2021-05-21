function [imgov,cmap,cscale,ovmax] = imgoverlay(img,overlay,ovmin,ovmax)
%
% [imgov cmap cscale ovmax] = imgoverlay(img,overlay,ovmin,<ovmax>)
%
% Creates an overlay given the underlay image (img),
% the overlay image (overlay) and the overlay display range
% (ovmin, ovmax).  This can be used to overlay statistics over
% a structural.  The img is rescaled to the range of 1:48, and
% the overlay image is rescaled to the range of 49:64, afterwhich
% the two are added together to create one image returned in 
% imgov.  Set the  colormap to cmap to display the underlay in
% gray scale and the overlay in color.  The overlay image can
% contain multiple planes.
%
%


%
% imgoverlay.m
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

if(nargin ~= 3 & nargin ~= 4)
  msg = 'USAGE: [imgov cmap cscale] = imgoverlay(img,overlay,ovmin,<ovmax>)'
  qoe(msg);error(msg);
end

if(~exist('imresize')) 
  fprintf(1,'INFO: Cannot find imresize function ...\n');
  fprintf(1,'      ... using fmri_imresize\n');
end

nImgRows  = size(img,1);
nImgCols  = size(img,2);
nImgDepth = size(img,3);

nOvRows  = size(overlay,1);
nOvCols  = size(overlay,2);
nOvDepth = size(overlay,3);

if(nImgDepth ~= 1 & nImgDepth ~= nOvDepth)
  msg = sprintf('Overlay (%d) and Underlay (%d) depths differ\n',...
                nImgDepth,nOvDepth);
  qoe(msg);error(msg);
end

% Rescale the image to the range 1:48 %
img = reshape1d(img);
img_min = min(img);
img_max = max(img);
img = ceil( 47 * (img - img_min)/ (img_max - img_min)) + 1;

overlay = reshape1d(overlay);

% Find max in overlay if it is not given %
if(nargin == 3)
  ovmax = max(overlay);
  iover = [];
else
  %% Find voxels above the max %%
  iover = find(overlay > ovmax);
end

%% Find voxels below the min %%
iunder = find(overlay < ovmin);

%% Clip the SupraMaxs to the Max %%
overlay(iover) = ovmax*ones(length(iover),1);

overlay = reshape(overlay, [nOvRows nOvCols nOvDepth]);
img     = reshape(img, [nImgRows nImgCols nImgDepth]);

for d = 1 : nOvDepth,

  if(nImgDepth > 1) k = d;
  else              k = 1;
  end

  if(~exist('imresize') | 1)
    tmp = reshape1d(fmri_imresize(overlay(:,:,d),[nImgRows nImgCols]));
  else
    tmp = reshape1d(imresize(overlay(:,:,d),[nImgRows nImgCols],'bicubic'));
  end

  iundermin = find(tmp<ovmin);
  iovermin  = find(tmp>ovmin);

  %% Compress the range of the overlays %%
  tmp = ceil( 15 * (tmp - ovmin)/ (ovmax-ovmin)) + 49;

  imgtmp = squeeze(img(:,:,k));
  tmp(iundermin) = imgtmp(iundermin);
  imgov(:,:,d) = reshape(tmp,[nImgRows nImgCols]);

end

% Construct the Colormap %
cmap = zeros(64,3);
cmap(1:48,1:3) = [[1:48]' [1:48]' [1:48]']/48;
j = jet(64);
cmap(49:64,:) = j(1:4:64,:);

% Construct the Color Scale %
dov = (ovmax-ovmin)/15;
cscale = [ovmin:dov:ovmax]';

return;
