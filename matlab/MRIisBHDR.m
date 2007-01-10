function [isBHDR, fstem] = MRIisBHDR(fspec)
% [isBHDR fstem] = MRIisBHDR(fspec)
%
% Determines whether the given file spec is an BHDR
% file based soley on its extension. BHDR is a bshort or
% bfloat.
%
% Returns non-zero if fspec is the name of an BHDR file,
%   Returns 1 if it has a .bhdr extension.
% Returns 0 otherwise.
%


%
% MRIisBHDR.m
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


fstem = [];
isBHDR = [];
if(nargin ~= 1)
  fprintf('[isBHDR fstem] = MRIisBHDR(fspec)\n');
  return;
end

isBHDR = 0;
if(length(fspec) < 6) return; end

ext = fspec(end-4:end);

if(strcmp(ext,'.bhdr')) 
  isBHDR = 1; 
  fstem = fspec(1:end-5);
end

return;













