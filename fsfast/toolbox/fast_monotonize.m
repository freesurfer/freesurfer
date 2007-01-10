function y = fast_monotonize(x)
% y = fast_monotonize(x)
%
% Makes columns of x monotonically decreasing
%
%


%
% fast_monotonize.m
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

if(nargin ~= 1)
  fprintf('y = fast_monotonize(x)\n');
  return;
end

nf = size(x,1);

y = x;
for n = 2:nf
  ind = find(y(n,:) > y(n-1,:));
  y(n,ind) = y(n-1,ind) ;
end

return;
