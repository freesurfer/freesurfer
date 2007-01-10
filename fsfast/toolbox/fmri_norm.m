function u = fmri_norm(v,order)
% u = fmri_norm(v, <order>)
% 
% normalizes the columns of v
%
% '$Id: fmri_norm.m,v 1.2 2007/01/10 22:02:33 nicks Exp $'


%
% fmri_norm.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
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

if(nargin == 0)
  msg = 'USAGE: u = fmri_norm(v, <order>)';
  qoe(msg);error(msg);
end

if(size(v,1)==1)
  u = ones(size(v));
  return;
end

if(nargin == 1 ) order = 1 ; end

f = (sum(v.^order)).^(1/order);
ind = find(f==0);
f(ind) = 10^10;
u = v ./ repmat(f,[size(v,1) 1]);

return;
