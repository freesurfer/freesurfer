function [tpr, tnc] = fast_glmpower(beta,X,rvar,C,alpha,nsides)
% [tpr tnc] = fast_glmpower(beta,X,rvar,C,alpha,nsides)
% 
% Computes the power (ie, true positive rate) of a hypothetical
% GLM analysis. NOTE: requires statistics toolbox.
%
% beta - vector of hypothetical regression coefficients
% X - deisgn matrix (dof = size(X,1) - size(X,2))
% rvar - hypothetical variance (scalar)
% C - contrast matrix (only one row allowed)
% alpha - Type I Error Rate (can be a vector)
% nsides - 1 for one-sided t test, or 2 for two-sided t test
%   default is one-sided.
%
% tpr = 1 - Type II Error Rate (vector same length as alpha)
% tnc = noncentrality parameter
%
% See also: fast_glmfitw, fast_fratiow, FTest.
% 
% $Id: fast_glmpower.m,v 1.1 2005/01/06 00:32:07 greve Exp $

tpr = [];
tnc = [];

if(nargin < 5 | nargin > 6)
  fprintf('[tpr tnc] = fast_glmpower(beta,X,rvar,C,alpha,nsides)\n');
  return;
end

if(size(C,1) ~= 1)
  fprintf('ERROR: C can only have 1 row\n');
  return;
end

if(~exist('nsides','var')) nsides = []; end
if(isempty(nsides)) nsides = 1; end

dof = size(X,1) - size(X,2);

% This is the noncentraliity parameter for the noncentral t distribution
tnc = C*beta/sqrt(rvar*C*inv(X'*X)*C');

% Compute the ideal TPR
if(nsides == 1)
  t = abs(tinv(alpha,dof)); % t corresonding to alpha in central t
else
  t = abs(tinv(alpha/2,dof)); % t corresonding to alpha in central t
end
tpr = 1 - nctcdf(t,dof,tnc); % get p from noncentral t

return;


