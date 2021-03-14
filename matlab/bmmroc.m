function [fpr,tpr,auc,fdr,nthresh,kappa] = bmmroc(pA,pI,lambda,ntrials,nthresh)
% [fpr tpr auc fdr nthresh kappa] = bmmroc(pA,pI,lambda,ntrials,<nthresh>)
%
% Computes the ROC (ie, TPR and FPR) based on a binomial mixture
% model. 
%
% pA - probability that a truly active voxel is declared positive (TPR)
% pI - probability that a truly inactive voxel is declared positive (FPR)
% lambda - mixure of active to total number of voxels
% ntrials 
% nthresh - threshold used to compute positives. Can be a vector. 
%  default is [0:ntrials].
%
% auc - area under the [ROC] curv
% fdr - false dicovery rate 
% kappa - Cohen's kappa (Thirion, NI, 2007)
%


%
% bmmroc.m
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


if(nargin < 4 | nargin > 5)
  fprintf('[fpr tpr auc fdr nthresh kappa] = bmmroc(pA,pI,lambda,ntrials,<nthresh>)\n');
  return;
end

if(~exist('nthresh','var')) nthresh = []; end
if(isempty(nthresh)) nthresh = [0:ntrials]; end

fpr = (1-binomialcdf(nthresh,ntrials,pI));
tpr = (1-binomialcdf(nthresh,ntrials,pA));

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
%den = (1-lambda)*fpr+lambda*tpr; % total number of positives
%ind = find(den == 0);
%den(ind) = 1;
%fdr = (1-lambda)*fpr./den;
%ind = find(fdr<0 | fdr>1);
%fdr(ind) = 0;
fdr = (1-lambda)*pI/(lambda*pA+(1-lambda)*pI);

% From Thirion, et al, 2007
piA1 = pA;    % fraction of active declared active
piA0 = 1-pA;  % fraction of active declared inactive
piI1 = pI;    % fraction of inactive declared active
piI0 = 1-pI;  % fraction of inactive declared inactive

% Fraction correctly classified as either active or inactive
p0 = lambda*piA1 + (1-lambda)*piI0;  
% Fraction declared inactive regardless of region. This is interpreted
% as the simple chance that a voxel is declared inactive, and 1-pi0 is
% the simple chance that a voxel is declared active.
pi0 = lambda*piA0 + (1-lambda)*piI0;
% Fraction of truly active voxels declared inactive by chance (chance errors)
p0A = lambda*pi0;
% Fraction of truly inactive voxels declared active by chance (chance errors)
p1I = (1-lambda)*(1-pi0);
% Fraction wrongly/incorrectly classified by chance
pW = p0A + p1I;  % pC in Thirion
% True (non-Thirion) Fraction correctly classified by chance
pCtrue = lambda*(1-pi0) + (1-lambda)*pi0;
% Correct classification rate minus chance wrong classification rate
kappaNum = (p0-pW);
% Best correct classification rate accounting for chance
kappaDen = (1-pW);
% Kappa
kappa = kappaNum/kappaDen;

return;








