function [zfil, inv_irf_fft] = fast_conv_invirf(z,irf)
% [zfil, inv_irf_fft] = fast_conv_invirf(z,irf)
% 
% Convolve input data with the inverse of the given impulse
% response function. z is nf-by-nc input data, where the filter
% will be applied across the row dimension (ie, nf). If the 
% length of irf is greather than nf, irf will be truncated.
% If the length of irf is less than nf, it will be padded with
% zeros. If greater than nf, it will be truncated. Uses FFT, 
% which is much faster than doing the convolution explicitly. 
% zfil will be the same size as z.
%


%
% fast_conv_invirf.m
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

zfil = [];
inv_irf_fft = [];

if(nargin ~= 2)
  fprintf('[zfil, inv_irf_fft] = fast_conv_invirf(z,irf)\n');
  return;
end

[nf nc] = size(z);

irf = reshape1d(irf);
if(length(irf) > nf) 
  % truncate down to nf%
  irf = irf(1:nf); 
elseif(length(irf) < nf) 
  % pad with zeros to reach nf
  irf = [irf; zeros(nf-length(irf),1)]; 
end

% Pad with zeros to avoid wrap-around %
irf = [irf; zeros(size(irf))];
z   = [z; zeros(size(z))];

% Compute inverse fft of irf %
irf_fft = fft(irf);
irf_fft_mag = abs(irf_fft);
irf_fft_ang = angle(irf_fft);
inv_irf_fft = exp(-irf_fft_ang)./(irf_fft_mag);

zfil_fft = fft(z) .* repmat(inv_irf_fft,[1 nc]);
zfil = real(ifft(zfil_fft));
zfil = zfil(1:nf,:);

return;
