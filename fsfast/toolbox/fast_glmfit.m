function [beta, rvar, dof] = fast_glmfit(y,X)
% [beta, rvar, dof] = fast_glmfit(y,X)
%

if(nargin ~= 2)
  fprintf('[beta, rvar, dof] = fast_glmfit(y,X)\n');
  return;
end

[nf nv] = size(y);
if(size(X,1) ~= nf)
  fprintf('ERROR: X and y have different number of frames\n');
  return;
end
dof = size(X,1) - size(X,2);

beta = (inv(X'*X)*X')*y;

if(nargout == 1) return; end

r = y - X*beta;
rvar = sum(r.^2)/dof;

return;



