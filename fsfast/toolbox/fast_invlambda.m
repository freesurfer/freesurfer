function invLambda = fast_invlambda(ErrCovMtx,nEigen)
% invLambda = fast_invlambda(ErrCovMtx,nEigen)


%
% fast_invlambda.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
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

if(nargin ~= 2)
  msg = 'USAGE: invLambda = fast_invlambda(ErrCovMtx,nEigen)';
  qoe(msg);error(msg);
end

NormErrCovMtx = fast_normcovmtx(ErrCovMtx);

[u s v] = svd(NormErrCovMtx);

nn = [1:nEigen];
ds = diag(s);
ds2 = ds(nn);
ids = 1./ds2;
is = diag(ids);
invLambda = u(:,nn) * is * u(:,nn)'; %'

return;
