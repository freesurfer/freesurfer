function z = fast_p2z(p)
% z = fast_p2z(p)
%
% Converts a significance to a z-score. p must be
% between -1 and 1
%
%
%


%
% fast_p2z.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/04/02 23:06:44 $
%    $Revision: 1.3 $
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

z = [];
if(nargin ~= 1) 
  fprintf('z = fast_p2z(p)\n');
  return;
end

z = sqrt(2)*erfcinv(2*abs(p)) .* sign(p);


return;










