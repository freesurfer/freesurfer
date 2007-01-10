function p = permutation(n,x,logflag)
% p = permutation(n,x,<logflag>);
% 
% Permutation of n and x:
% p = (n) =    n!
%     (x)   x!(n-x)!
%
% n is an integer > 0
% x is an integer > 0 and < n
%
% If logflag = 1, then returns log(p) instead of p
%
% Uses logarithms to avoid overflow.
%
%


%
% permutation.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
%    $Revision: 1.4 $
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

p = [];
if(nargin ~= 2 & nargin ~= 3)
  fprintf('p = permutation(n,x,<logflag>)\n');
  return;
end
if(~exist('logflag','var')) logflag = 0; end

if(round(n) ~= n)
  fprintf('ERROR: n = %g, an integer\n',n);
  return;
end

if(round(x) ~= x)
  fprintf('ERROR: x = %g, an integer\n',x);
  return;
end

if(n < 0)
  fprintf('ERROR: n = %d, must be >= 0\n',n);
  return;
end

if(x < 0 | x > n)
  fprintf('ERROR: x = %d, must be >= 0 and <= n\n',x);
  return;
end

% Use logarithms to avoid overflows %
lp = sum(log(1:n)) - sum(log(1:x)) - sum(log(1:n-x));
if(logflag) p = lp;
else        p = exp(lp);
end

return;






