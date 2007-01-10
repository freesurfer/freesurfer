function slice = fast_ldslice(volid,sliceno)
% slice = fast_ldslice(volid,sliceno)


%
% fast_ldslice.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
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

slice = [];

if(nargin ~= 2)
  msg = 'USAGE: slice = fast_ldslice(volid,sliceno)'
  qoe(msg) ; error(msg);
end

fmt = fast_getvolformat(voldid);
if(isempty(fmt))
  msg = sprintf('Could not determine format of %s',volid);
  qoe(msg) ; error(msg);
end

switch(fmt)

end



return;
