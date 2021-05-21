function f = pulsespectrum(x0,x1,w)
% f = pulsespectrum(x0,x1,w)
%
% Continuous (analytical) spectrum of a pulse starting 
% at x0 and ending at x1. w is the radial frequency.
%
% If x0 and x1 have more than one element, then a spectrum
% is produced, and f will have multiple columns.
% 
% If x0 = -x1, then the result will be all real.
%
% Example:
%   ph = 2*pi * [0:N-1]'/N; ph = ph - ph(N/2+1); 
%     ph will go from -pi to (almost) pi, passing through 0
%     at N/2+1. This makes the x axis go from -N/2 to (N-1)/2
%     passing through 0 at N/2 + 1: x = [-N/2 : (N-1)/2]; So let
%   x0 = -10; x1 = 10; then
%   f = pulsespectrum(x0,x1,ph);
%   p = fftshift(abs(ifft(fftshift(f))));
%   plot(x,p) shows the pulse. The pulse starts at about -10 and
%   ends at about +10. There is a lot of ringing and roll off
%   (ie, it's not a crisp pulse). I assume that this is because 
%   the spectrum vector is sampled from the continuous, which is 
%   not what the FFT assumes. Incidentally, the decoding DFT matrix
%   is iDFT = exp(+i*(ph*x))/N; ie, p = abs(iDFT*f); This p will
%   be the same as the one above. Less ringing when x0,x1 are on
%   the half voxel (they are continuous).
%
% abs(f) =  sqrt(2*(1-cos(w.*(x1-x2))))./w;
%
% d(abs(f))/dw = -a./(w.^2) + dx*sin(dx*w)./(w.*a);
%  a = sqrt(2*(1-cos(dx*w)));
%
%


%
% pulsespectrum.m
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

f = [];

if(nargin ~= 3)
  fprintf('f = pulsespectrum(x0,x1,w)\n');
  return;
end

if(length(x0) ~= length(x1))
  fprintf('ERROR: x0 and x1 have different sizes\n');
  return;
end

d = x1-x0;
ind = find(d <= 0);
if(~isempty(ind))
  fprintf('ERROR: x0 >= x1\n');
  return;
end
nx = length(x1);

indwz  = find(w==0);
indwnz = find(w~=0);

f = zeros(prod(size(w)),nx);
for nthx = 1:nx
  x00 = x0(nthx);
  x11 = x1(nthx);

  fx = zeros(size(w));
  fx(indwz) = x11-x00;

  wx0 = w(indwnz)*x00;
  wx1 = w(indwnz)*x11;
  fx(indwnz) = (cos(wx1) - cos(wx0) + i*(sin(wx0)-sin(wx1)) )./(-i*w(indwnz));

  f(:,nthx) = reshape1d(fx);

end

f = reshape(f,[size(w) nx]);

return
