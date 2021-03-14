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
%
%


%
% fast_polytrendmtx.m
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
  r = r0 - Xptm*(inv(Xptm'*Xptm)*Xptm'*r0);
  r = r/std(r);
  Xptm = [Xptm r];
end

Xptm_run = zeros(ntrs,nruns*(order+1));
n1 = (run-1)*(order+1) + 1;
n2 = n1 + order;
Xptm_run(:,n1:n2) = Xptm;


if(nargout == 2)
  Rptm = eye(ntrs) - Xptm*inv(Xptm'*Xptm)*Xptm';
end

return;
