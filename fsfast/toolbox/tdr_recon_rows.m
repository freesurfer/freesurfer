function [rimg, dgbeta, ncsub] = tdr_recon_rows(kimg,Rrow,PERev)
% [rimg, dgbeta] = tdr_recon_rows(kimg,Rrow,<PERev>)
%
% The rows are reconned by applying Rrow. Includes deghosting.
%
% See also tdr_recon_cols.
% 
% $Id: tdr_recon_rows.m,v 1.4 2003/12/19 22:25:23 greve Exp $

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
    [rimg(:,:,slice,frame) dgbeta(:,slice,frame)] = ...
	tdr_deghost(kslice,Rrow,PERev);
  end
end

return;












