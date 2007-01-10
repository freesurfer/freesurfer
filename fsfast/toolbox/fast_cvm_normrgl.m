function [mrgl, niters, alpha] = fast_cvm_normrgl(m,condmin)
% mrgl = fast_cvm_normrgl(m,condmin)
%
% Regularization and normalization.
%


%
% fast_cvm_normrgl.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
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

nitersmax = 1000;
mrgl = [];
if(nargin ~= 1 & nargin ~= 2)
  msg = 'USAGE: mrgl = fast_cvm_normrgl(m,condmin)';
  qoe(msg);error(msg);
end

if(nargin == 1)
  condmin = inf;
end

if(condmin < 1) 
  msg = 'Cannot specify a condition less than 1';
  qoe(msg);error(msg);
end

mrgl = m;
mineig = min(eig(mrgl));
mcond = cond(mrgl);
ddm = diag(diag(m));
alpha = .1;
niters = 1;
while(niters < nitersmax & (mineig < 0 | mcond > condmin))
  mrgl = fast_cvm_normalize(m + alpha*ddm);
  mineig = min(eig(mrgl));
  mcond = cond(mrgl);
  alpha = alpha + .1;
  fprintf('%4d  %6.4f   %g  %g\n',niters,alpha,mineig,mcond);
  niters = niters + 1;
end

return;
