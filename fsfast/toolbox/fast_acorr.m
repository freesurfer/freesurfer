function R = fast_acorr(x,scaling,dof)
% R = fast_acorr(x,<scaling>,<dof>)
%
% x is the input data (nf-by-nc)
% scaling is the scaling option
%   'none' - returns sum of the cross products at each delay
%   'biased' - divides none by nf
%   'unbiased' - divides none by nf-lag
%   'unbiasedcoeff' - divides unbiased by zero-lag (careful!)
%   'coeff' - divides biased by zero-lag <default>
% dof is the number of degrees of freedom in x. Default is nf.
%   dof only has an effect when unibasedcoeff is chosen. dof can 
%   be handy when computing the ACF of residuals. Note: the ACF
%   at lags >= dof is invalid.
%
% Computes the autocorrelation of x using an FFT.  This should
% produce the same results as the native matlab function xcorr
% when invoked as xcorr(x,size(x,1),scaling).
%
% There are a few differences between fast_acorr and xcorr:
%   1. fast_acorr does not compute cross-correlations
%   2. fast_acorr will compute a separate ACF for each column in x
%   3. the ACF from fast_acorr will be one-sided, ie, R(1) is for 
%      zero delay, R(2) is for delay 1, etc. 
%   4. fast_acorr has an unbiasedcoeff option. Note: this option
%      can produce acfs that are greater than 1 (eventhough
%      it is supposed to be a coefficient).
%   5. fast_acorr can incorporate dof information
%   6. fast_acorr is much faster because it uses an FFT.
% 
% $Id: fast_acorr.m,v 1.3 2003/04/15 03:51:32 greve Exp $

if(nargin ~= 1 & nargin ~= 2 & nargin ~= 3)
  msg = 'USAGE: R = fast_acorr(x,<scaling>)';
  qoe(msg);error(msg);
end

% Get dims of x %
[ntrs ncols] = size(x);

if(~exist('scaling')) scaling = []; end
%if(isempty(scaling))  scaling = 'unbiasedcoeff'; end
if(isempty(scaling))  scaling = 'coeff'; end
if(~exist('dof'))     dof = ntrs; end


% FFT of X %
fftx = fft(x,2*ntrs);

% Correlation in time is element-by-conj(element) mult in frequency %
% Note: Convolution would not use conj().
fftx2 = fftx .* conj(fftx); 
clear fftx;

% Inverse FFT to get back to time-domaine
R = real(ifft(fftx2));
clear fftx2;

% Only keep the first ntrs
R = R(1:ntrs,:);

% At this point R has no scaling

if(strcmp(scaling,'unbiased') | strcmp(scaling,'unbiasedcoeff'))
  lag = [0:ntrs-1]'; %'
  c = dof - lag;
  icz = find(c==0);
  c(icz) = 1;
  R = R./repmat(c,[1 ncols]);
end

if(strcmp(scaling,'biased'))
  R = R/ntrs;
end

if(strcmp(scaling,'coeff') | strcmp(scaling,'unbiasedcoeff'))
  % Divide by the zeroth lag
  R0 = R(1,:);
  iz = find(R0==0);
  R0(iz) = 10^10;
  R = R./repmat(R0,[ntrs 1]);
end

return;
