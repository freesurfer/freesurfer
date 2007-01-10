function siteno = fbirn_siteno(sitename)
% siteno = fbirn_siteno(sitename)
%
% Returns the site number id given the site name
% 
%


%
% fbirn_siteno.m
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

if(nargin ~= 1)
  fprintf('siteno = fbirn_siteno(sitename)\n');
  return;
end

siteno = [];

sitelist = fbirn_sitelist;
nsites = size(sitelist,1);
for nthsite = 1:nsites
  if(strcmp(sitename,deblank(sitelist(nthsite,:))))
    siteno = nthsite;
    return;
  end
end


return;  
  






