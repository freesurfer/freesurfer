function [item, count] = sscanfitem(tline,nthitem)
% [item count] = sscanfitem(tline,nthitem)
%
% Reads the nth item from a string list of items separated by white
% space. The item is returned as a string. If successful, count=1,
% otherwise count=0 (and item is empty).
%
%


%
% sscanfitem.m
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

item = '';
count = 0;

if(nargin ~= 2)
  fprintf('[item count] = sscanfitem(tline,nthitem)\n');
  return;
end

fmt = '%s';
for n = 1:nthitem-1
  fmt = sprintf('%%*s %s',fmt);
end

[item count] = sscanf(tline,fmt,1);

return;





