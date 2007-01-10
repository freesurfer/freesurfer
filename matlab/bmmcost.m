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


%
% bmmcost.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%


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

