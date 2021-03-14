function [avgRxx, stdRxx] = fmri_acorr(x, nMaxLag, ind) 
%
% [avgRxx stdRxx] = fmri_acorr(x)
% [avgRxx stdRxx] = fmri_acorr(x, nMaxLag) 
% [avgRxx stdRxx] = fmri_acorr(x, nMaxLag, ind) 
%
% Computes the autocorrelation of x using an FFT.  This should
% produce the same results as the native matlab function xcorr
% when invoked as xcorr(x,nMaxLag,'unbiased').  This function
% should be much faster because it uses an FFT.  Rxx is scaled so 
% that the zeroth lag is 1. If x is a matrix, the Rxx of each column 
% is computed separately, and the average and std at each delay point
% are returned.
%
% The input argument ind contains a list of columns in x to use
% when computing the autocorrelation.
%
% Size of x: nRows, nCols, nTP, nRuns
%
%


%
% fmri_acorr.m
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


if(nargin ~= 1 & nargin ~= 2 & nargin ~= 3)
  msg = 'USAGE: [avgRxx stdRxx] = fmri_acorr(x,<<nMaxLag>,ind>)';
  qoe(msg);error(msg);
end

% Get dims of x %
[nRows nCols Ntp nRuns] = size(x);
nVoxels = nRows*nCols;

x = permute(x, [3 1 2 4]);
x = reshape(x, [Ntp nVoxels nRuns]);

% Set the maximum lag depending upon whether nMaxLag is specified %
if(nargin == 2 | nargin == 3 )   nML = nMaxLag;
else                             nML = Ntp-1;
end

% Check that the max lag is not too large
if(nML > Ntp-1)
  msg = 'nMaxLag must be <= Ntp in waveform';
  qoe(msg);error(msg);
end

for r = 1:nRuns,

  % autocorrelation using fft %
  if(nargin == 3)
    fftVEW = fft(x(:,ind(:,r),r),2*Ntp);
  else
    fftVEW = fft(x(:,:,r),2*Ntp);
  end

  fftVEW2 = fftVEW .* conj(fftVEW);

  %% matrix of autocor funtions , Ntp x Nv %%
  Rsum = real(ifft(fftVEW2)); 
  Rsum = Rsum(1:nML+1,:); % extract only the components of interest.

  % DOF at each lag
  dof = Ntp - [0:nML]'; %'
  dof2 = repmat(dof, [1 size(Rsum,2)]);

  % correct for dof at each lag (unbiased)
  Rvox = Rsum./dof2; 

  % Scale for unity at zero lag %
  maxRvox = repmat(max(Rvox), [nML+1 1]);
  Rvox = Rvox./maxRvox;

  % Average across voxels %
  avgR1xx = mean(Rvox,2);
  stdR1xx = std(Rvox,[],2);

  % Make two-sided %
  avgRxx(:,r) = [avgR1xx(nML+1:-1:2)' avgR1xx']'; %'
  stdRxx(:,r) = [stdR1xx(nML+1:-1:2)' stdR1xx']';

end

return;
