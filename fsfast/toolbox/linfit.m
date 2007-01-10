function [slope, offset] = linfit(y,x)
%
% [slope offset] = linfit(y,x)
%
% Solves the equation:
%   y = slope * x + offset
%
%


%
% linfit.m
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

if(nargin ~= 2)
  msg = 'USAGE: [slope offset] = linfit(y,x)';
  qoe(msg);error(msg);
end 

if(length(y) ~= length(x))
  msg = 'Length of x and y must be the same';
  qoe(msg);error(msg);
end

x = reshape1d(x);
y = reshape1d(y);

x2 = [x ones(size(x))];

m = pinv(x2)*y;

slope  = m(1);
offset = m(2);

return;
