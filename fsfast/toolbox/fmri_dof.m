function dof = fmri_dof(X,DTOrder)
%
% Compute degrees of freedom.
%
% dof = fmri_dof(X,DTOrder)
%
% X - stim conv matricies: nTP x nTotEst x nRuns, where nTotEst
%     is (nHEst + nPreStim)*nNonNullCond.
% DTOrder - number of trends removed.
%
% The dof is the number of observations (nRuns*nTP) minus
% the number of parameters fit (HOrder*nNNCond + nRuns*PPOrder).
%
% $Id: fmri_dof.m,v 1.1 2003/03/04 20:47:39 greve Exp $
%

if(nargin ~= 2)
  msg = 'USAGE: dof = fmri_dof(X,DTOrder)';
  qoe(msg);
  error(msg);
end

[nTP nTotEst nRuns] = size(X);

dof = nRuns*nTP - nRuns*DTOrder - nTotEst;

return;
