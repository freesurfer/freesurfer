function r = rande(rdim,order)
%
% r = rande(dim,<order>)
%
% Generates matrix of given dimension whose elements 
% have an erlang distribution. The expected value
% is 1. The variance is 1/r. The stddev is 1/sqrt(r).
%
% If the order is unspecified, defalts to 1.
% An order of 1 is an exponential distribution.
%
% r = rande([5 50],3);
%
% See also erlang.m
%
%


%
% rande.m
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

if(exist('rdim')  ~= 1) rdim  = 1; end
if(exist('order') ~= 1) order = 1; end

r = 0;
for n = 1:order
  r = r + -log(rand(prod(rdim),1));
end
r = r/order; % make avg=1
r = reshape(r, [rdim 1]);

return;
