function [yest] = fmri_estsignal(X,hest)
%
% yest = fmri_estsignal(X,hest)
%
% X:     Stim Conv Mtx (nTP x nTotEst x nRuns)
% hest:  hemodynamic resp est  (nRows x nCols x nTotEst)
% 
% yest:  estimated fMRI slice (nRows x nCols x nTP x nRuns)
%
% $Id: fmri_estsignal.m,v 1.1 2003/03/04 20:47:39 greve Exp $
%

if(nargin ~= 2)
  msg = 'USAGE: [yest] = fmri_estsignal(X,hest)';
  qoe(msg);
  error(msg);
end

nTP     = size(X,1);
nTotEst = size(X,2);
nRuns   = size(X,3);

nRows   = size(hest,1);
nCols   = size(hest,2);

if(size(hest,3) ~= nTotEst)
  msg = 'hest and X dimensions are inconsistent';
  qoe(msg);
  error(msg);
end

% Compute number of voxels %
nV = nRows*nCols; 

hest = reshape(hest, [nV nTotEst])';

for r = 1:nRuns,
  yest(:,:,r) = X(:,:,r) * hest;
end

yest = permute(yest, [2 1 3]);
yest = reshape(yest, [nRows nCols nTP nRuns]);

return;
