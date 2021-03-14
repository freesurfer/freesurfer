function R = fast_acorr(x,scaling,dof,tpexclude)
% R = fast_acorr(x,<scaling>,<dof>,<tpexclude>)
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
% tpexclude - 1-based time points to exclude from the calculation. 
%   The values of x at those time points are set to 0, and
%   the scaling is adjusted at each delay to account for the
%   fact that not as many data points go into the calculation.
%   This scaling occurs regardless of the scaling argument.
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
%


%
% fast_acorr.m
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

R = [];
if(nargin < 1 | nargin > 4)
  msg = 'USAGE: R = fast_acorr(x,<scaling>,<dof>,<tpexclude>))';
  qoe(msg);error(msg);
end

% Get dims of x %
[ntrs ncols] = size(x);

if(~exist('scaling','var'))   scaling = []; end
if(isempty(scaling))          scaling = 'coeff'; end
if(~exist('dof','var'))       dof = ntrs; end
if(~exist('tpexclude','var')) tpexclude = []; end

% Set points to exclude equal to 0 (more tpexc stuff below)
if(~isempty(tpexclude))
  % But first make sure they are in range
  if(max(tpexclude) > ntrs)
    fprintf('ERROR: tpexclude max (%d) > ntrs (%d)\n', ...
	    max(tpexclude),ntrs);
    return;
  end
  if(min(tpexclude) < 1)
    fprintf('ERROR: tpexclude min (%d) < 1\n',min(tpexclude));
    return;
  end
  x(tpexclude,:) = 0;  
end

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

% Apply scaling for tpexclude
if(~isempty(tpexclude))
  tpall = ones(ntrs,1);
  tpinclude = tpall;
  tpinclude(tpexclude) = 0;
  sall  = fast_acorr(tpall,'none');
  sincl = fast_acorr(tpinclude,'none');
  indz  = find(sincl==0);
  %fprintf('nz = %d\n',length(indz));
  indnz = find(sincl~=0);
  sc = ones(ntrs,1);
  sc(indnz) = sall(indnz)./sincl(indnz);
  R = R ./repmat(sc,[1 ncols]);
  if(~isempty(indz))  R(indz,:) = 0; end
end

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
