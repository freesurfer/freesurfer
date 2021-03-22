function [beta, rvar, t, tsig] = fast_vvglm(y,x,polyorder)
% [beta, rvar, t, tsig] = fast_vvglm(y,x,<polyorder>)
%
% Performs a voxel-by-voxel GLM analysis in which each column of x
% and y are treated as their own set of equations to solve, ie,
%   y(:,c) = x(:,c)*beta(c) + n
%
%


%
% fast_vvglm.m
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

beta = [];

if(nargin < 2 | nargin > 3)
  fprintf('[beta rvar t tsig] = fast_vvglm(y,x,<polyorder>)\n');
  return;
end

if(~exist('polyorder','var')) polyorder = []; end

[nfy nvy] = size(y);
[nfx nvx] = size(y);

if(nfy ~= nfx)
  fprintf('ERROR: x and y do not have same number of frames\n');
  return;
end

if(nvy ~= nvx)
  fprintf('ERROR: x and y do not have same number of voxels\n');
  return;
end

nv = nvx;
nf = nfx;
dof = nf - 1;

% Detrend both, does not bias beta values
if(~isempty(polyorder))
  Xdt = fast_polytrendmtx(1,nf,1,polyorder);
  Rdt = eye(nf) - Xdt*inv(Xdt'*Xdt)*Xdt';
  x = Rdt*x;
  y = Rdt*y;
  dof = dof - polyorder;
end


% Now do voxel-by-voxel fit
sumx2 = sum(x.^2);
indz = find(sumx2<eps);
sumx2(indz) = 1e10;
sumxy = sum(x.*y);
beta = sumxy ./ sumx2;
yhat = x .* repmat(beta,[nf 1]);
r = y - yhat;
rvar = sum(r.^2)/dof;
indz = find(rvar<eps);
rvar(indz) = 1e10;

% Do t-test
betavar = rvar ./ sumx2;
t = beta./sqrt(betavar);
tsig = tTest(dof,t);


% This is to test linear FPRs under null hypothesis
% [pdf, alpha, nxhist, fpr] = ComputePDF(tsig,.01,1,.01);
% loglog(alpha,alpha,alpha,fpr)

return;




% This was just to verify that this gives the same results as my
% other programs. It does.
ytmp = y(:,1);
Xtmp = x(:,1);
[betatmp rvartmp vdof] = fast_glmfit(ytmp,Xtmp);
[F Fsig] = fast_fratio(betatmp,Xtmp,rvartmp,1);

keyboard

return;





