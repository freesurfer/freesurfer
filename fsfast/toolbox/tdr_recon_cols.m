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
%    $Date: 2007/01/10 22:02:35 $
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

img = [];

if(nargin ~= 1)
  fprintf('img = tdr_recon_rows(kimg)\n');
  return;
end

%[nrows ncols nslices nframes] = size(kimg);
img = fftshift(ifft(fftshift(kimg,1),[],1),1);

return;

