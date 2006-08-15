function [fpr,tpr,auc] = bmmroc(pA,pI,lambda,ntrials,thresh)
% [fpr tpr auc] = bmmroc(pA,pI,lambda,ntrials,thresh);
%
% Computes the ROC (ie, TPR and FPR) based on a binomial mixture
% model. 
%
% $Id: bmmroc.m,v 1.1 2006/08/15 23:42:51 greve Exp $

Nthresh = round(ntrials*thresh);
fpr = 1-binomialcdf(Nthresh,ntrials,pI);
tpr = 1-binomialcdf(Nthresh,ntrials,pA);
%auc = trapz(fpr,tpr);
auc = 0;

return;








