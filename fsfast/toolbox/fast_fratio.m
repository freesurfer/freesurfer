function [F, dof1, dof2, g] = fast_fratio(beta,X,rvar,C)
% [F, dof1, dof2, g] = fast_fratio(beta,X,rvar,C)
%
% p = FTest(dof1, dof2, F, <dof2max>)
%

if(nargin ~= 4)
  fprintf('[F, dof1, dof2, g] = fast_fratio(beta,X,rvar,C)\n');
  return;
end

gcvm = inv(C*inv(X'*X)*C');
J = size(C,1);
nv = size(beta,2);
dof1 = J;
dof2 = size(X,1)-size(X,2);

indz  = find(rvar == 0);
indnz = find(rvar ~= 0);
if(~isempty(indz))
  beta = beta(:,indnz);
  rvar = rvar(:,indnz);
end

g = C*beta;
F = (sum(g .* (gcvm*g))./rvar)/J;

if(~isempty(indz))
  F0 = zeros(1,nv);
  F0(indnz) = F;
  F = F0;
end

return;



