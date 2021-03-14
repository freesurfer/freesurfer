function [d1, d2, d3, match, vxc] = fast_esttrans(v1,v2,frame)
% [d1, d2, d3, match, vxc] = fast_esttrans(v1,v2,frame)
%
% Estimate translation between two volumes. The result is 
% the amount that v2 must be translated to best match v1.
% match is a number between 0 and 1 indicating how good
% the match is.
% 
% Example:
% n = 64;
% v1 = randn(n,n,n);
% v2 = fast_mshift(v1, [ 1 2 3 ], 1);
% [d1 d2 d3] = fast_esttrans(v1,v2);
% v1b = fast_mshift(v2, [ d1 d2 d3 ], 1);
% v1b will now be idendical to v1


%
% fast_esttrans.m
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

frame = 1;

v1fft = fftn(v1(:,:,:,frame));
v2fft = fftn(v2(:,:,:,frame));

vxcfft = v1fft .* conj(v2fft);
vxc = abs(ifftn(vxcfft));

% Rescale so that the peak is between 0 and 1 %
tmp = (sum(reshape1d(v1.^2)) + sum(reshape1d(v2.^2)))/2;
vxc = vxc/tmp;

[match ivxcmax] = maxmd(vxc);

d1 = ivxcmax(1)-1;
d2 = ivxcmax(2)-1;
d3 = ivxcmax(3)-1;

if(d1 > size(v1,1)/2) d1 =  d1 - size(v1,1) ; end
if(d2 > size(v1,2)/2) d2 =  d2 - size(v1,2) ; end
if(d3 > size(v1,3)/2) d3 =  d3 - size(v1,3) ; end


return;
