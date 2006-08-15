function cost = bmmcost(params,Ndet)
% cost = bmmcost(params,Ndet)
%  
% cost (log-likelihood) function for binomial mixture model
%
% params = [pA pI lambda]
% cost = -llbmm(Ndet,pA,pI,lambda)
%
% optparams = fminsearch('bmmcost',initparams,[],Ndet);
%
% Init is important because there is a symmetry between the
% parameters, ie, you get the same cost if you swap pA and pI and use
% 1-lambda.
%
% $Id: bmmcost.m,v 1.1 2006/08/15 23:41:28 greve Exp $

pA = params(1);
pI = params(2);
lambda = params(3);

% Check for out-of-bounds
if(pA < eps | (1-pA) < eps | ...
   pI < eps | (1-pI) < eps | ...
   lambda < eps | (1-lambda) < eps)
  cost = 10e10;
  return;
end

% Use negative because it is a minimization
cost = -llbmm(Ndet,pA,pI,lambda);

return

