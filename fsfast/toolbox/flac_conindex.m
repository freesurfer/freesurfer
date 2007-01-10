function ind = flac_conindex(conname,flac)
% ind = flac_conindex(conname,flac)
% Returns the index of the given contrast name in the flac
%


%
% flac_conindex.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
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

ind = [];
if(nargin ~= 2)
  fprintf('ind = flac_conindex(conname,flac)\n');
  return;
end

ncon = length(flac.con);
for nthcon = 1:ncon
  if(strcmp(flac.con(nthcon).name,conname))
    ind = nthcon;
    return;
  end
end

% Will return empty ind if gets here.

return;
















