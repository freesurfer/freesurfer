function [K, phmat] = fast_dftmtx(N_or_Phase)
% [K phmat] = fast_dftmtx(N_or_Phase)
% 
% Computes the matrix that implements an N-point 1D DFT.
%
% N_or_Phase can be one of two things:
%   N - the number of points in the DFT (scalar)
%   Phase - the phase associated with each sample (vector)
%
% Specifying N is equivalent to specifying Phase = 2*pi*[0:N-1]/N';
%
% For a vector f = randn(N,1); then F = fftshift(K*f); will be the
% same as fft(f). Actually, only the magnitudes will be the same.
%
% $Id: fast_dftmtx.m,v 1.2 2003/09/10 03:00:56 greve Exp $

if(nargin ~= 1)
  fprintf('[K phmat] = fast_dftmtx(N_or_Phase)\n');
  return;
end

if(length(N_or_Phase)==1) 
  N = N_or_Phase;
  ph0 = 2*pi * [0:N-1]'/N;
else                      
  N = length(N_or_Phase);
  ph0 = reshape(N_or_Phase,[N 1]);
end

r = [0:N-1];
r = r - r(N/2 + 1);
ph0 = ph0 - ph0(N/2+1);
phmat = (ph0 * r); % outer product
K = exp(-i * phmat);

return;