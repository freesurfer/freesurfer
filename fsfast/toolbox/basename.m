function name = basename(path)
% name = basename(path)
% 
% This is an attempt to recreate the unix basename 
% command in matlab.
%


%
% basename.m
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

name = [];

if(nargin ~= 1)
  msg = 'USAGE: name = basename(path)'
  qoe(msg); error(msg);
end

len = length(path);

if(len == 1) 
  if(strcmp(path,'/'))
    name = '/';
  else
    name = path;
  end
  return;
end

% strip trailing '/' character %
if(strcmp(path(len),'/'))
  path = path(1:len-1);
  len = len-1;
end

for n = len-1:-1:1
  if(strcmp(path(n),'/'))
    name = path(n+1:len);    
    return;
  end
end

name = path;

return;
