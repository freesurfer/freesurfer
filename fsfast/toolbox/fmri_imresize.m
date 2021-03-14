function b = fmri_imresize(a, bdim, method)
% b = fmri_imresize(a, bdim)
% b = fmri_imresize(a, bdim, method)
%
% Resize input image a to bdim.
%
% a is the input image
% b is the output image
% bdim is the dimensions (rows,cols) of the output
%   which must be an integral multiple of the input dim.
% method - interpoltation method (currently only 'nearest'
%   is suppored (this is the default).
%
%


%
% fmri_imresize.m
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

if(nargin ~= 2 & nargin ~= 3)
  msg = 'USAGE: b = fmri_imresize(a, bdim, <method>)';
  qoe(msg);error(msg);
end

[arows acols aimgs] = size(a);

brows = bdim(1);
bcols = bdim(2);

frow = brows/arows;
if(frow ~= round(frow))
  msg = sprintf(...
   'Output rows (%d) is not an integer mult of Input  rows (%d)',...
   brows,arows);
  qoe(msg);error(msg);
end

fcol = bcols/acols;
if(fcol ~= round(fcol))
  msg = sprintf(...
   'Output cols (%d) is not an integer mult of Input  cols (%d)',...
   bcols,acols);
  qoe(msg);error(msg);
end

avoxs = arows*acols;
a0 = a;
a = reshape(a, [1 avoxs*aimgs]);

b1 = repmat(a,   [frow 1]);
b2 = reshape(b1, [brows acols aimgs]);
b2 = permute(b2, [2 1 3]);

b3 = reshape(b2, [1 brows*acols*aimgs]);
b4 = repmat(b3,   [fcol 1]);
b5 = reshape(b4, [bcols brows aimgs]);
b5 = permute(b5, [2 1 3]);
b  = b5;

if(0)
figure(1);
subplot(3,2,1);
imagesc(a0);
subplot(3,2,2);
imagesc(b2);
subplot(3,2,3);
imagesc(b3);
subplot(3,2,4);
imagesc(b4);
subplot(3,2,5);
imagesc(b5);
end

return
