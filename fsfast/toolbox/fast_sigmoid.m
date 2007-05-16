function y = fast_sigmoid(x,dx,s);
% y = fast_sigmoid(x,<dx>,<s>);
%
% y passes thru 0.5 at x=dx. Default dx=0.
% s controls the slope. Default s=1.
% dx and s can either be scalars or same size as x.
%
% y = 0 at x = -inf 
% y = 1 at x = +inf 
% y will be the same size as x. 
%
% $Id: fast_sigmoid.m,v 1.1 2007/05/16 06:08:06 greve Exp $

%
% fast_sigmoid.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/05/16 06:08:06 $
%    $Revision: 1.1 $
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

y = [];
if(nargin < 0 | nargin > 3)
  fprintf('y = fast_sigmoid(x,<dx>,<s>);\n');
  return;
end

if(~exist('dx','var')) dx = []; end
if(isempty(dx)) dx = 0; end

if(~exist('s','var')) s = []; end
if(isempty(s)) s = 1; end

y = 1./(1 + exp(-s.*(x-dx)));

return;
