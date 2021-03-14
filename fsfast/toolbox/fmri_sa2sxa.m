function [hAvg, eVar, hStd] = fmri_sa2sxa(ySA,nHEst)
%
% [hAvg eVar hStd] = fmri_sa2sxa(ySA,nHEst)
%
% Converts from selavg format to selxavg format. In selavg
% format, the averages and standard deviations are interleaved
% in the same data structure for all conditions, including
% fixation.  
%
% -------------- Input Arguments ---------------
% 1. ySA - interleaved averages and standard deviations for
%      all conditions (incl fix) and all delays.
%      size(ySA) = [nRows nCols 2*nDelays*nCond].
% 2. nHEst - number of coefficients to estimate in the hemodynamic
%      response.
%
% ------------- Output Arguments ---------------
% 1. hAvg - the average hemodynamic impulse response(s) for
%      all delays and non-null stimulus conditions
%      (size(hAvg) = [nRows nCols nHEst*nNNCond]).
% 2. eVar - the varaince of the residual errors on a voxel by
%      voxel basis (size(eVar) = [nRows nCols]).
% 1. hStd - the standard deviations of the  hemodynamic response(s) for
%      all delays and non-null stimulus conditions
%      size(hStd) = [nRows nCols nHEst*nNNCond]).
%
%


%
% fmri_sa2sxa.m
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

if(nargin ~= 2 & nargin ~= 3)
  msg = 'USAGE: [hAvg, eVar, hStd] = fmri_sa2sxa(ySA,nHEst)';
  qoe(msg);
  error(msg);
end

% Extact relevant dimensions %
nRows   = size(ySA,1);
nCols   = size(ySA,2);
Nv      = nRows*nCols;
Nas     = size(ySA,3); % includes fixation
nCond   = Nas/(2*nHEst); % includes fixation
Nnnc    = nCond -1; % number of non-null conditions
Nch     = Nnnc*nHEst;

eVar = (ySA(:,:,nHEst+1).^2); % grab a single plane %

ySA = permute(ySA,[3 1 2]);
ySA = reshape(ySA,[nHEst 2 nCond nRows nCols]);

hAvg = squeeze(ySA(:,1,[2:nCond],:,:));
hAvg = reshape(hAvg, [Nch nRows nCols]);
hAvg = permute(hAvg, [2 3 1]);

hStd = squeeze(ySA(:,2,[2:nCond],:,:));
hStd = reshape(hStd, [Nch nRows nCols]);
hStd = permute(hStd, [2 3 1]);

return;

%%%%%%%%%% Ignore everything below this line %%%%%%%%%%%%%%%%


%% Get Cond 0 to Subtract from the others %%%
%% This has no effect on selxavg format but
%% is necessary for selavg format.
hAvg0 = squeeze(ySA(:,1,1,:,:));
hAvg0 = permute(hAvg0, [3 4 1 2]);
hAvg0 = reshape(hAvg0,[nRows nCols nHEst]);
hAvg0 = repmat(hAvg0,[1 1 Nnnc]);

hAvg = squeeze(ySA(:,1,[2:nCond],:,:));
hAvg = permute(hAvg, [3 4 1 2]);
hAvg = reshape(hAvg,[nRows nCols Nch]);
hAvg = hAvg - hAvg0;

hStd = squeeze(ySA(:,2,[2:nCond],:,:));
hStd = permute(hStd, [3 4 1 2]);
hStd = reshape(hStd,[nRows nCols Nch]);

return;

%%%%%%%%%%%% Ignore stuff below this line %%%%%%%%%%%%


ind = [1:nHEst];
i = [];
for n = 2:nCond,
  i = [i (ind+2*nHEst*(n-1))];
end

hAvg = ySA(:,:,i);

if(nargin ~= 3)
  DOF = ones(nCond,1)/nCond; % assume all equal %
end
DOF = reshape(DOF, [1 1 nCond])./ DOF(1);

mDOF = zeros(nRows,nCols,(nCond-1)*nHEst);
for n = 1:nCond-1,
  %fprintf('%d %d %d\n',n,nHEst*(n-1)+1,nHEst*n-1);
  mDOF(:,:,[nHEst*(n-1)+1:nHEst*n]) = repmat(DOF(n), [nRows nCols nHEst]);
end
mDOF0 = DOF(1) * ones(size(mDOF));

% Subtract off the zero condition.  For selxavg, this has
% no effect.  For selavg, this makes the hAvg comensurate
% with selxavg.
hAvg0 = repmat(ySA(:,:,[1:nHEst]),[1 1 size(hAvg,3)/nHEst]);
hAvg = hAvg - hAvg0;

if(nargout == 3)
  hStd  = ySA(:,:,i+nHEst);
  hStd0 = ySA(:,:,[nHEst+1:2*nHEst]);
  hStd0 = repmat(hStd0,[1 1 size(hAvg,3)/nHEst]);
  hVar  = hStd.^2;
  hVar0 = hStd0.^2;
  hStd  = sqrt(hVar .* mDOF + hVar0 .* mDOF0);
end


return;
