function [fpr,tpr,auc,fdr,nthresh] = bmmroc(pA,pI,lambda,ntrials,nthresh)
% [fpr tpr auc fdr nthresh] = bmmroc(pA,pI,lambda,ntrials,<nthresh>)
%
% Computes the ROC (ie, TPR and FPR) based on a binomial mixture
% model. 
%
% pA - probability that a truly active voxel is declared positive
% pI - probability that a truly inactive voxel is declared positive
% lambda - mixure of active to total number of voxels
% ntrials 
% nthresh - threshold used to compute positives. Can be a vector. 
%  default is [0:ntrials].
%
% auc - area under the [ROC] curv
% fdr - false dicovery rate as a function of threshold
%
% $Id: bmmroc.m,v 1.2 2006/08/16 03:17:13 greve Exp $

if(nargin < 4 | nargin > 5)
  fprintf('[fpr tpr auc] = bmmroc(pA,pI,lambda,ntrials,<nthresh>)\n');
  return;
end

if(~exist('nthresh','var')) nthresh = []; end
if(isempty(nthresh)) nthresh = [0:ntrials]; end

fpr = 1-binomialcdf(nthresh,ntrials,pI);
tpr = 1-binomialcdf(nthresh,ntrials,pA);

% Area under the ROC curve
auc = abs(trapz(fpr,tpr)); % needs abs for reversal

% False dicovery rate
den = ((1-lambda)*fpr+lambda*tpr);
ind = find(den == 0);
den(ind) = 1;
fdr = fpr./den;
ind = find(fdr<0 | fdr>1);
fdr(ind) = 0;

return;








