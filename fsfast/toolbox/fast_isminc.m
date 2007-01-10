function isbvol = fast_isbvol(volid)
% isbvol = fast_isbvol(volid)
%
% Returns 1 if volume is in bfile format, ie, there exist file
% with name volid_%03d.bshort or .bfloat.  Actually, it just
% looks for the header volid_%03d.hdr.
% 
%


%
% fast_isminc.m
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

isbvol = 0;

if(nargin ~= 1)
  msg = 'USAGE: r = fast_isbvol(volid)'
  qoe(msg); error(msg);
end

stem = deblank(volid);

nslices = 0;
for slice = 0:30
  fname = sprintf('%s_%03d.hdr',stem,nslices);
  fid = fopen(fname,'r');
  if(fid ~= -1) 
    fclose(fid);
    isbvol = 1;
    return;
  end
end

return;
