function r = fast_fileexists(filename)
% r = fast_fileexists(filename)
% 1 if it exists and is readable , 0 if not


%
% fast_fileexists.m
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

if(nargin ~= 1) 
  msg = 'USAGE: r = fast_fileexists(filename)';
  qoe(msg); error(msg);
end

fid = fopen(filename,'r');
if(fid == -1 ) r = 0;
else
  r = 1;
  fclose(fid);
end

return;
