function sitename = fbirn_sitename(siteno)
% sitename = fbirn_sitename(siteno)
%
% Returns the site name given the site number.
% 
%


%
% fbirn_sitename.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
%    $Revision: 1.3 $
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

sitename = [];

if(nargin ~= 1)
  fprintf('sitename = fbirn_sitename(siteno)\n');
  return;
end

sitelist = fbirn_sitelist;
nsites = size(sitelist,1);

if(siteno > nsites)
  fprintf('ERROR: siteno = %d > nsites = %d\n',siteno,nsites);
  return;
end

sitename = deblank(sitelist(siteno,:));

return;  
  






