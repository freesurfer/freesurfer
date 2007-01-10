function r = direxists(dirpath)
% r = direxists(dirpath)
%
% returns 1 if dirpath exists, 0 else
%
%


%
% direxists.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:29 $
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

r = -1;
if(nargin ~= 1)
  fprintf('r = direxists(dirpath)\n');
  return;
end

d = dir(dirpath);
if(size(d,1) == 0) r = 0;
else               r = 1;
end

return;





