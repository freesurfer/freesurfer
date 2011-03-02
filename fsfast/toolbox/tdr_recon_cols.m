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
%


%
% tdr_recon_cols.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:07 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

img = [];

if(nargin ~= 1)
  fprintf('img = tdr_recon_rows(kimg)\n');
  return;
end

%[nrows ncols nslices nframes] = size(kimg);
img = fftshift(ifft(fftshift(kimg,1),[],1),1);

return;

