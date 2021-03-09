function [rimg, dgbeta, ncsub] = tdr_recon_rows(kimg,Rrow,PERev)
% [rimg, dgbeta] = tdr_recon_rows(kimg,Rrow,<PERev>)
%
% The rows are reconned by applying Rrow. Includes deghosting.  kimg
% can be any number of slices and frames. Flips even rows prior to
% recon.
%
% See also tdr_recon_cols.
% 
%


%
% tdr_recon_rows.m
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

rimg = [];
dgbeta = [];

if(nargin ~= 2 & nargin ~= 3)
  fprintf('[rimg, dgbeta] = tdr_recon_rows(kimg,Rrow,<PERev>)\n');
  return;
end
if(exist('PERev') ~= 1) PERev = 0; end

[nrows ncols nslices nframes] = size(kimg);
evenrows = [2:2:nrows];

rimg = zeros(nrows,size(Rrow,2),nslices,nframes);
dgbeta = zeros(2,nslices,nframes);
ncsub = zeros(nslices,nframes);
for frame = 1:nframes
  for slice = 1:nslices
    kslice = kimg(:,:,slice,frame);
    kslice(evenrows,:) = fliplr(kslice(evenrows,:));
    %if(perev) kepi = flipud(kepi);  end 
    [rimg(:,:,slice,frame) dgbeta(:,slice,frame)] = ...
	tdr_deghost(kslice,Rrow,PERev);
  end
end

return;












