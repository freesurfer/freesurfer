function img = tdr_recon_cols(kimg)
% img = tdr_recon_rows(kimg)
%
% The columns are reconned by applying the inverse FFT.
%
% kimg is the complex k-space image (may have multiple slices,
% frames). Most of the time, the rows of kimg will have already 
% been reconned with tdr_recon_rows in order to remove ghosting.
%
% The abs() is not taken.
%
% See also tdr_recon_rows, tdr_kshift.
% 
% $Id: tdr_recon_cols.m,v 1.1 2003/10/28 04:33:13 greve Exp $

img = [];

if(nargin ~= 1)
  fprintf('img = tdr_recon_rows(kimg)\n');
  return;
end

%[nrows ncols nslices nframes] = size(kimg);
img = fftshift(ifft(fftshift(kimg,1),[],1),1);

return;

