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
% Example: if v is a stack of images in the order of acquisition,
%  then so = tdr_sliceorder(nslices,2) and v2 = v(:,:,so), then
%  v2 will be the same stack but with adjacent anatomical slices
%  next to each other.
% 
% $Id: tdr_sliceorder.m,v 1.1 2003/10/20 22:15:53 greve Exp $

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
  so = [1:2:nslices 2:2:nslices];
  return;
end

M = floor(nslices/2) + 1;
c = [1:M; M+1:2*M];
so = reshape1d(c);
so = so(1:nslices);

return;










