function [y, H] = randnorm(m, n, Sigma)
%
% [y H] = randnorm(m, n, Sigma)
%
% Creates an m-by-n matrix y whose columns are
% distributed N(0,Sigma).  H is the cholesky 
% decomposition of Sigma. y = H * randn(m,n).
%
%


%
% randnorm.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
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

if(nargin ~= 3) 
  msg = 'USAGE: [y H] = randnorm(m, n, Sigma)';
  qoe(msg); error(msg);
end

q = size(Sigma,1) ; 
if( q > m)
  tmp = Sigma(1:m,1:m);
  Sigma = tmp;
  clear tmp;
elseif(q < m)
  tmp = zeros(m,m);
  tmp(1:q,1:q) = Sigma;
  Sigma = tmp;
  clear tmp;
end

H = chol(Sigma);

y = H * randn(m,n);

return ;
