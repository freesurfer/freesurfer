function vclip = fast_clip(v,thresh,clipdir)
% vclip = fast_clip(v,<<thresh>,clipdir>)
%
% clipdir: <positive>, negative, both
%
%


%
% fast_clip.m
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

if(nargin < 1  | nargin > 3) 
  msg = 'vclip = fast_clip(v,<<thresh>,clipdir>)';
  qoe(msg);error(msg);
end

if(~exist('thresh') | isempty(thresh))  thresh  = 0; end
%if(isempty(thresh)) thresh  = 0; end

if( ~exist('clipdir') | isempty(clipdir) ) clipdir = 'positive'; end

if(~strncmp(clipdir,'both',1) & ~strncmp(clipdir,'positive',1) & ...
   ~strncmp(clipdir,'negative',1))
  msg = sprintf('clipdir = %s, must be either pos, neg, or both',clipdir);
  qoe(msg);error(msg);
end

vclip = v;

if(strncmp(clipdir,'both',1) | strncmp(clipdir,'positive',1))
  ind = find(v>thresh);
  vclip(ind) = thresh;
end

if(strncmp(clipdir,'both',1) | strncmp(clipdir,'negative',1))
  ind = find(v < -thresh);
  vclip(ind) = -thresh;
end

return;
