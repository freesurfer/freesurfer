function volformat = fast_getvolformat(volid)
% volformat = fast_getvolformat(volid)
% 
%


%
% fast_getvolformat.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
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

volformat = [];

if(nargin ~= 1)
  msg = 'USAGE: volformat = fast_getvolformat(volid)'
  qoe(msg); error(msg);
end

if(fast_ismincvol(volid))  volformat = 'minc';
elseif(fast_isbvol(volid)) volformat = 'bfile';
end


return;
