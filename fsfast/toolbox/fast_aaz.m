function [aaz, fstd] = fast_aaz(fslice,dtorder,fftflag)
% [aaz fstd] = fast_aaz(fslice,<dtorder>,<fftflag>)
%
% Computes the absolute z score averaged over spatial voxels.
%
% fslice is the functional data and can be either
%  nrows-ncols-nframes or nframes-nvoxels. The FFT
%  can only be done on the latter.
% dtorder - temporal detrending order. Default is 0
%  (ie, remove the mean). Use -1 for no detrending.
% fftflag - apply a 2D sptial fft to each slice.
%  Default is no fft.
% aaz - Average Absolute Z, one for each frame.
%
%
%
% (c) Douglas N. Greve, 2004
%


%
% fast_aaz.m
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

aaz = [];

if(nargin < 1 | nargin > 3)
  fprintf('[aaz fstd] = fast_aaz(fslice,dtorder,fftflag)\n');
  return;
end

if(~exist('dtorder','var')) dtorder = []; end
if(isempty(dtorder)) dtorder = 0; end

if(~exist('fftflag','var')) fftflag = []; end
if(isempty(fftflag)) fftflag = 0; end

szfslice = size(fslice);

% Do a spatial 2D FFT
if(fftflag) 
  if(length(szfslice) ~= 3)
    fprintf('ERROR: fslice must be 2D + Time with FFT flag\n');
    return;
  end
  fslice = abs(fft2(fslice)); 
end

if(length(szfslice) == 3)
  % Reshape if necessary
  % 2D + Time
  [nr nc nf] = size(fslice);
  nv = nr*nc;
  fslice = reshape(fslice,[nv nf])';
else
  % Time by Spatial Voxels
  [nf nv] = size(fslice);
end

% Detrend %
if(dtorder > -1)
  X = fast_polytrendmtx(1,nf,1,dtorder);
  beta = (inv(X'*X)*X')*fslice;
  fslice = fslice - X*beta;
end

% Compute functional stddev
fstd = std(fslice);

% Check for nothing there, return null
indnz = find(fstd > eps);
nnz = length(indnz);
if(nnz == 0) return; end 

% Remove voxels with zero stddev
if(nnz ~= nv)
  fslice = fslice(:,indnz);
  fstd = fstd(indnz);
end

% Compute Z by dividing by stddev
z = fslice ./ repmat(fstd,[nf 1]);

% Compute Average Abs Z, average over space
aaz = mean(abs(z),2);

return;
