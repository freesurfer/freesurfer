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


%
% bmmroc.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.5 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


if(nargin < 4 | nargin > 5)
  fprintf('[fpr tpr auc] = bmmroc(pA,pI,lambda,ntrials,<nthresh>)\n');
  return;
end

if(~exist('nthresh','var')) nthresh = []; end
if(isempty(nthresh)) nthresh = [0:ntrials]; end

fpr = 1-binomialcdf(nthresh,ntrials,pI);
tpr = 1-binomialcdf(nthresh,ntrials,pA);

% This is a bit of a hack to account for the fact that the cdf
% gives the probability of nhits below or *at* a given threshold
% (since the threshold is discrete). So when the threshold is 0,
% the probability is non-zero which means that the fpr is not
% 1. This adjusts it.
fpr(2:end) = fpr(1:end-1);
fpr(1) = 1;
tpr(2:end) = tpr(1:end-1);
tpr(1) = 1;

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








