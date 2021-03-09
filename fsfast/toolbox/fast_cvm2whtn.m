function [W, r, s, nunder ] = fast_cvm2whtn(cvm,nmax,pctrmagmin)
% [W, r, s, nunder] = fast_cvm2whtn(cvm,nmax,pctrmagmin)
%
% Computes whitening matrix from a covariance matrix. 
%
% Algorithm: computes the autocorrelation function (r) as the
% average of the diagonal components of the CVM. This 
% toeplitzization (sorry) of the CVM enforces the assumption
% that the underlying process is time-invariant. The magnitude
% and phase of the FFT of the autocorrelation function is computed.
% The magnitude of the whitening filter is then computed as the
% inverse of the square root of that of the autocorrelation 
% function. The phase of the whitening filter is computed as the
% negative of that of the autocorrelation function. The whitening
% filter is then converted back to the time domain to become the
% whitening function (s). This is then used to compute an upper 
% triangular toeplitz matrix, the whitening matrix (W), which has 
% the same size as cvm.
%
% If nmax is specified, then only the first nmax components of
% the whitening function are used (the rest are padded to zero).
%
% If pctrmagmin is specified, then the the frequency components of
% the magnitude of the FFT of the autocorrelation function are not
% allowed to go below Rmagmin = max(Rmag)*pctrmagmin/100, ie, a
% a percentage of maximum magnitude. The number that fall under this
% threshold is returned as nunder.
%
%


%
% fast_cvm2whtn.m
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

% Testing: 
% ntrs = 100; ncols = 1000; alpha = .5; rho = .7;
% nmax = 20; pctmin = 10;
% y = synthcornoise(ntrs,ncols,alpha,rho);
% ycvm = y*y'; %'
% [W r s nunder] = fast_cvm2whtn(ycvm,nmax,pctmin);  
% ycvm2 = W*ycvm*W'; %'
% [W2 r2 s2 nunder2] = fast_cvm2whtn(ycvm2);
% plot(2:nmax,r(2:nmax),2:nmax,r2(2:nmax));

W = [];
if(nargin < 1 | nargin > 3)
  msg = 'USAGE: [W, r, s] = fast_cvm2whtn(cvm,nmax,pctrmagmin)';
  qoe(msg);error(msg);
end

ntrs = size(cvm,1);

if(nargin < 2 | isempty(nmax))       nmax = ntrs; end
if(nargin < 3 | isempty(pctrmagmin)) pctrmagmin = 0; end

if(nmax > ntrs)
  msg = sprintf('nmax (%d) > ntrs (%d)\n',nmax,ntrs);
  qoe(msg);error(msg);
end

if(pctrmagmin < 0 | pctrmagmin >= 100)
  msg = sprintf('pctrmagmin = %g, must be bet 0 and 100\n',pctrmagmin);
  qoe(msg);error(msg);
end

% Estimate autocorrelation function as the mean of the various
% diagonal elements.
r = zeros(ntrs,1);
for n = 1:ntrs
  r(n) = mean(diag(cvm,n-1));
end

% Normalize by the first component %
r = r/r(1);

% Convert autocorrelation function to frequency %
R      = fft(r(1:nmax));
Rmag   = abs(R);
Rphase = angle(R);

% Make sure that Rmag does not go below a minimum %
Rmagmax = max(Rmag);
Rmagmin = Rmagmax*pctrmagmin/100;
iunder = find(Rmag < Rmagmin);
nunder = length(iunder);
Rmag(iunder) = Rmagmin;

% Compute the frequency mag and phase of the whitening filter
Smag = 1./sqrt(Rmag);
Sphase = -Rphase;

% Compute real/imaginary components of whitening filter
S = Smag .* cos(Sphase) + sqrt(-1)*Smag .* sin(Sphase);

% Convert back to time
s0 = real(ifft(S));

% Normalize so that the first component is always 1
s0 = s0/s0(1);

% Pad with zeros to get back to ntrs %
s = zeros(ntrs,1);
s(1:nmax) = s0(1:nmax);

% Construct upper triangular toeplitz matirx 
tmp = zeros(ntrs,1);
tmp(1) = 1;
W = toeplitz(tmp,s);

return;
