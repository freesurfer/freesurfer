function ulist = stringunique(strvlist)
% ulist = stringunique(strvlist)
%
% Returns a unique list of strings found in the vertically
% concatenated string list strvlist. Matches are determined
% AFTER deblanking.
%
% $Id: stringunique.m,v 1.1 2007/10/08 23:46:45 greve Exp $

%
% stringunique
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/10/08 23:46:45 $
%    $Revision: 1.1 $
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


ulist = [];
if(nargin ~= 1)
  fprintf('ulist = stringunique(strvlist)\n');
  return;
end

nlist = size(strvlist,1);

% Go thru each string in the list
for nth = 1:nlist
  nthstr = deblank(strvlist(nth,:));
  hit = 0;
  % Go thru each string in the unique list
  for mth = 1:size(ulist,1)
    ustr = deblank(ulist(mth,:));
    if(strcmp(nthstr,ustr))
      % The nth string is already in the unique list
      hit = 1;
      break;
    end
  end
  if(hit == 0)
    % The nth string was not found in the unique list, so add it
    ulist = strvcat(ulist,nthstr);
  end
end

return;
