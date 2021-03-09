function [F, Fsig, ces, cesvar, edof] = fast_fratio(beta,X,rvar,C,Sn,dof2max,edofuse)
% [F, Fsig, ces, cesvar] = fast_fratio(beta,X,rvar,C,<Sn>,<dof2max>,<edofuse>)
%
% Sn is the covariance matrix of the noise AFTER any filtering.
% edofuse is the DOF to use instead of that computed form the
%  various inputs. This may be needed if projections were done
%  on the raw time course data prior to analysis.
%
% Worsley, K.J. and Friston, K.J. Analysis of fMRI Time-Series
% Revisited - Again. Neuroimage 2, 173-181, 1995.
%
% See also: fast_glmfit.m
%
%


%
% fast_fratio.m
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

if(nargin < 4 | nargin > 7)
  fprintf('[F, Fsig, ces, cesvar, edof] = fast_fratio(beta,X,rvar,C,<Sn>,<dof2max>,<edofuse>)\n');
  return;
end

nf = size(X,1);

if(exist('dof2max') ~= 1) dof2max = []; end
if(exist('Sn') ~= 1)      Sn = []; end
if(exist('edofuse') ~= 1) edofuse = []; end

if(~isempty(Sn))
  if(size(Sn,1) ~= nf)
    fprintf('ERROR: Sn dimension mismatch\n');
    return;
  end
  R = eye(nf)-X*inv(X'*X)*X';
  vdof = trace(R*Sn);
  edof = (vdof.^2)/trace(R*Sn*R*Sn);
  % Covariance matrix of contrast effect size
  cescvm = inv(C*inv(X'*X)*X'*Sn*X*inv(X'*X)*C');
else
  if(isempty(edofuse))
    edof = size(X,1) - size(X,2);
  else
    edof = edofuse;
  end
	
  % Covariance matrix of contrast effect size
  cescvm = inv(C*inv(X'*X)*C');
end

J = size(C,1);
nv = size(beta,2);
dof1 = J;
dof2 = edof;

if(J==1) cesvar = rvar./cescvm; end

if(nv > 1)
  indz  = find(rvar == 0);
  indnz = find(rvar ~= 0);
  if(~isempty(indz))
    beta = beta(:,indnz);
    rvar = rvar(:,indnz);
  end
else
  indz = [];
  indnz = 1;
end


% Contast Effect Size
ces = C*beta;

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
    ces0 = zeros(J,nv);
    ces0(:,indnz) = ces;
    ces = ces0;
    clear Fsig0 ces0;
  end
end

return;



