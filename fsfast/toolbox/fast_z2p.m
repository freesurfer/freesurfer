function p = fast_z2p(z)
% p = fast_z2p(z)
% converts z values into p values (ie, the area under the curve 
% to the right of the z value). This is a one-tailed test (ie,
% the p-value does not take the sign of the z).
%
%


%
% fast_z2p.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
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

p = erfc((z(:))/sqrt(2))/2;
p = reshape(p,size(z));

% Do this if you are signing the z
%p = p.*sign(z);
%ind = find(z==0);
%p(ind) = 1;

return;




