function [Xptm_run, Rptm] = fast_polytrendmtx(run,ntrs,nruns,order)
% [Xptm Rptm] = fast_polytrendmtx(run,ntrs,nruns,order)
%
% Polynomial trend matrix. 
%  order=0 - mean offset
%  order=1 - mean + linear
% All columns are orthonormal. The first (mean) column is always 1
% so that it's regression coefficient is the mean offset.
%
% If nruns > 1, then X is padded horizontally with zeros to account
% for extra runs in the design matrix.
%
% If Rptm is specified, then the residual forming matrix for a
% single run is computed as:
%   Rptm = eye(ntrs) - X*inv(X'*X)*X';
%
% $Id: fast_polytrendmtx.m,v 1.3 2004/06/03 20:07:05 greve Exp $
%

if(nargin ~= 4)
  msg = 'USAGE: Xptm = fast_polytrendmtx(run,ntrs,nruns,order)';
  qoe(msg);error(msg);
end

if(order < 0) 
  Xptm_run = [];
  Rptm = [];
  return;
end

Xptm = ones(ntrs,1);
t = [0:ntrs-1]'; %'
for n = 1:order
  r0 = t.^n;
  M = eye(ntrs) - Xptm*inv(Xptm'*Xptm)*Xptm';
  r = M*r0;
  r = r/std(r);
  Xptm = [Xptm r];
end

Xptm_run = zeros(ntrs,nruns*(order+1));
n1 = (run-1)*(order+1) + 1;
n2 = n1 + order;
Xptm_run(:,n1:n2) = Xptm;

if(nargout)
  Rptm = eye(ntrs) - Xptm*inv(Xptm'*Xptm)*Xptm';
end

return;
