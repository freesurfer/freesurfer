function isMGH = MRIisMGH(fspec)
% isMGH = MRIisMGH(fspec)
%
% Determines whether the given file spec is an MGH/MGZ
% file based soley on its extension
%
% Returns non-zero if fspec is the name of an MGH file,
%   Returns 1 if it has a .mgh extension.
%   Returns 2 if it has a .mgz extension.
% Returns 0 otherwise.
%


%
% MRIisMGH.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
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


isMGH = [];
if(nargin ~= 1)
  fprintf('isMGH = MRIisMGH(fspec)\n');
  return;
end

isMGH = 0;
if(length(fspec) < 5) return; end

ext = fspec(end-3:end);

if(strcmp(ext,'.mgh')) isMGH = 1; end
if(strcmp(ext,'.mgz')) isMGH = 2; end

return;













