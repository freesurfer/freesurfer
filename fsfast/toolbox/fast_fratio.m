function [F, Fsig, ces, edof] = fast_fratio(beta,X,rvar,C,Sn,dof2max)
% [F, Fsig, ces] = fast_fratio(beta,X,rvar,C,<Sn>,<dof2max>)
%
% Sn is the covariance matrix of the noise after any filtering.

if(nargin < 4 | nargin > 6)
  fprintf('[F, Fsig, ces, edof] = fast_fratio(beta,X,rvar,C,<Sn>,<dof2max>)\n');
  return;
end

nf = size(X,1);

if(exist('dof2max') ~= 1) dof2max = []; end
if(exist('Sn') ~= 1)      Sn = []; end

if(~isempty(Sn))
  if(size(Sn,1) ~= nf)
    fprintf('ERROR: Sn dimension mismatch\n');
    return;
  end
  R = eye(nf)-X*inv(X'*X)*X';
  vdof = trace(R*Sn);
  edof = (vdof.^2)/trace(R*Sn*R*Sn);
else
  edof = size(X,1) - size(X,2);
end

J = size(C,1);
nv = size(beta,2);
dof1 = J;
dof2 = edof;

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



