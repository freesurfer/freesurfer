function [C, nC] = fmri_pcacovmtx(Q,p)
%
% [C, nC] = fmri_pcacovmtx(Q,p)
%
% Computes a new covariance matrix based on the principal components
% of the given matrix Q.  The number of components used will be
% enough to explain at least 100*p percent of the variance.
%
%
% $Id: fmri_pcacovmtx.m,v 1.1 2003/03/04 20:47:40 greve Exp $
%

n = size(Q,1);
nRuns = size(Q,3);
M = zeros(size(Q));

for r = 1:nRuns,

  %fprintf(1,'-- PCA for run %d\n',r);
  M = Q(:,:,r);

  [MEigVect MEigVal] = eig(M);
  EigVal = diag(MEigVal);

  % Fraction of variance explained by each eigenvector
  fExplain = EigVal/sum(EigVal); 

  % Sort fractions in ascending order.
  [fExpSort indSort] = sort(fExplain); 

  % Reverse so in descending order 
  fExpSort = fExpSort(n:-1:1);
  indSort  = indSort(n:-1:1);

  % Cumulative Fraction explained
  cfExpSort = cumsum(fExpSort);

  % Components needed to explain p
  nC(r) = min(find(cfExpSort>p));

  % Construct a matrix of these eigenvectors
  MEigVect2 = MEigVect(:,indSort(1:nC));

  % Construct a matrix of these eigenvalues
  EigVal2  = ones(n,1) * min(EigVal(indSort(1:nC)));
  EigVal2(indSort(1:nC)) = EigVal(indSort(1:nC));
  MEigVal2 = diag(EigVal2);

  % Recompute the matrix with subset of eigen vectors%
  C(:,:,r) = MEigVect * MEigVal2 * MEigVect';

end

return;
