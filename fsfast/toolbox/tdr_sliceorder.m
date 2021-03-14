function so = tdr_sliceorder(nslices,which)
% so = tdr_sliceorder(nslices,which)
%
% Slice order utility for interleaved slice aquisition
%
% which = 1, convert an aquired slice number into an anatomical
%  slice number. anatsliceno = so(acqsliceno)
%
% which = 2, convert an anat slice number into an acquired
%  slice number. acqsliceno = so(anatsliceno)
%
% Siemens changes the interleaved slice order depending upon
% whether there are an even or odd number of slices, eg
%   odd(5):  1 3 5 2 4
%   even(6): 2 4 6 1 3 5
%
% Example: if v is a stack of images in the order of acquisition,
%  then so = tdr_sliceorder(nslices,2) and v2 = v(:,:,so), then
%  v2 will be the same stack but with adjacent anatomical slices
%  next to each other.
% 
%


%
% tdr_sliceorder.m
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

so = [];
if(nargin ~= 2)
  fprintf('so = tdr_sliceorder(nslices,which)\n');
  return;
end

if(which ~= 1 & which ~= 2)
  fprintf('ERROR: which = %d, must be 1 or 2\n',which);
  return;
end

if(which == 1)
  % converts acq slice to anat slice
  if(rem(nslices,2)==1) % odd
    so = [1:2:nslices 2:2:nslices];
  else % even
    so = [2:2:nslices 1:2:nslices ];
  end
  return;
end

so0 = tdr_sliceorder(nslices,1);
[tmp so] = sort(so0);

return

%M = floor(nslices/2) + 1;
%c = [1:M; M+1:2*M];
%so = reshape1d(c);
%so = so(1:nslices);
%return;










