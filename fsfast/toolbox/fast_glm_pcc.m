function rho = fast_glm_pcc(beta,X,C,rvar)
% rho = fast_glm_pcc(beta,X,C,rvar)
%
% Computes the partial correlation coefficient of the given
% univariate contrast.
%
% WARNING: X must have a column of 1s or both X and y must have
% been demeaned before the glm.
%

if(nargin ~= 4) 
  fprintf('rho = fast_glm_pcc(beta,X,C,rvar)\n');
  return;
end

[J nbeta] = size(C);
if(J ~= 1)
  fprintf('ERROR: fast_glm_pcc: contrast must be univariate\n');
  return;
end

ntp = size(X,1);
DOF = ntp - nbeta;

% Design matrix projected onto the contrast space. This forms
% the space of the contrast
Xc = X*C';

% Null space of the contrast
D = null(C);
% Design matrix projected onto the null-space contrast. This forms
% the space of the "nuisance" regressors
Xd = X*D;
% Residual-forming matrix of the nuisance regressors
Rd = eye(ntp) - Xd*inv(Xd'*Xd)*Xd';
% Orthogonalize both Xc and yhat wrt the nuisance regressors
Xcd = Rd*Xc;
yhatd = Rd*(X*beta);

% Compute sums and cross-products needed for correlation coef
Xcdtyz = Xcd'*yhatd;
sumXcd = sum(Xcd);
sumyz = sum(yhatd); % Fails without column of 1s in X
XcdtXcd = sum(Xcd.^2);
yztyz = sum(yhatd.^2) + DOF*rvar;

rho = (Xcdtyz-sumXcd.*sumyz)./sqrt((XcdtXcd-sumXcd.^2).*(yztyz-sumyz.^2));

return;
