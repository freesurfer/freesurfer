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
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
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
