function [dsm, g] = fast_smooth1d(d,fwhm)
% [dsm g] = fast_smooth1d(d,<fwhm>)
%
% Gaussian smooths the columns of d. Pads with zeros so no
% wrap-around. Uses fft. Should work on complex data. 
%
%


%
% fast_smooth1d.m
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

% Test 1
% N = 50; gstd = 2;
% y = zeros(N,1); y(round(N/2)) = 1;
% [ysm g] = fast_smooth1d(y,gstd*sqrt(log(256.0)););
% x = [1:N]'; x = x - x(round(N/2));
% g0 = gaussian(x,0,gstd);
% plot(x,g0,'o-',x,ysm);
%
% Test 2
% nn = 1:N;
% y = randn(N,3000);
% ysm = fast_smooth1d(y,gstd);
% ymss = mean(ysm.^2,2);
% plot(nn,ymss,nn,g0'*g0)

dsm = [];
if(nargin < 1 | nargin > 2)
  fprintf('[dsm g] = fast_smooth1d(d,<fwhm>)\n');
  return;
end

if(~exist('fwhm','var')) fwhm = []; end
if(isempty(fwhm)) fwhm = 1; end
gsigma = fwhm/sqrt(log(256.0));

nrows = size(d,1);
ncols = size(d,2);
nrows2 = 2*nrows;

x = ([1:nrows2]' - (nrows2/2 + 1))/gsigma;
g = exp( -(x.^2)/2 );
g = g/sum(g);

dpad = zeros(nrows2,ncols);
dpad(1:nrows,:) = d;

kd = fft(dpad);
kg = fft(g);

dsm = ifft(kd .* repmat(kg,[1 ncols]));
dsm = fftshift(dsm,1); % so you can take 1:nrows
dsm = dsm(1:nrows,:);

if( max(abs(imag(d(:)))) == 0 )
  % Only take the real if input is totally real
  dsm = real(dsm);
end


return;
