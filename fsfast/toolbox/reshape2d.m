function m = reshape2d(x)
% m = reshape2d(x)
% Reshape into a 2d array where
% size(m,1) = size(x,1)


%
% reshape2d.m
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

if(nargin ~= 1)
  msg = 'm = reshape2d(x)'
  qoe(msg);error(msg);
end

szx = size(x);
xdim = length(szx);

if(xdim < 2)
  msg = sprintf('Dimension of x is only %d',xdim);
  qoe(msg);error(msg);
end

m = reshape(x, [szx(1) prod(szx(2:xdim))]);
return
