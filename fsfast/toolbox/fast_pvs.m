function [pvs, u, ds] = fast_pvs(y)
% [pvs, u, s] = fast_pvs(y)
% 
% Computes percent variance spanned by each of the eigenvectors
%


%
% fast_pvs.m
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

pvs = [];
u = [];
ds = [];

if(nargin ~= 1)
  fprintf('[pvs, u, s] = fast_pvs(y)\n');
  return;
end

My = y*y'; %'
[u s v] = svd(My);
ds = diag(s);
pvs = 100*ds/sum(ds);

return;
