function [F, Fsig, ces] = fast_fratio(beta,X,rvar,C,dof2max)
% [F, Fsig, ces] = fast_fratio(beta,X,rvar,C,<dof2max>)
%

if(nargin < 4 | nargin > 5)
  fprintf('[F, Fsig, ces] = fast_fratio(beta,X,rvar,C,<dof2max>)\n');
  return;
end
if(exist('dof2max') ~= 1) dof2max = []; end

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

% Contast Effect Size
ces = C*beta;

% Covariance matrix of contrast effect size
cescvm = inv(C*inv(X'*X)*C');

if(J ~= 1) F = (sum(ces .* (cescvm*ces))./rvar)/J;
else       F = ((ces.^2)./rvar)*(cescvm/J);
end

if(nargout > 1)
  Fsig = FTest(dof1, dof2, F, dof2max);
end

if(~isempty(indz))
  F0 = zeros(1,nv);
  F0(indnz) = F;
  F = F0;
  clear F0;
  if(nargout > 1)
    Fsig0 = ones(1,nv);
    Fsig0(indnz) = Fsig;
    Fsig = Fsig0;
    clear Fsig0;
  end
end

return;



