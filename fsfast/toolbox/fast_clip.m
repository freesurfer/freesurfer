function vclip = fast_clip(v,thresh,clipdir)
% vclip = fast_clip(v,<<thresh>,clipdir>)
%
% clipdir: <positive>, negative, both
%
% $Id: fast_clip.m,v 1.1 2003/03/04 20:47:37 greve Exp $

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
