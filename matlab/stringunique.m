function ulist = stringunique(strvlist)
% ulist = stringunique(strvlist)
%
% Returns a unique list of strings found in the vertically
% concatenated string list strvlist. Matches are determined
% AFTER deblanking.
%
% $Id: stringunique.m,v 1.2 2011/03/02 00:04:13 nicks Exp $

%
% stringunique
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:13 $
%    $Revision: 1.2 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
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
