function isminc = fast_ismincvol(volid)
% isminc = fast_ismincvol(volid)
%
% Returns 1 if volume is in MINC format.
% 
%


%
% fast_ismincvol.m
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

isminc = 0;

if(nargin ~= 1)
  msg = 'USAGE: r = fast_ismincvol(volid)'
  qoe(msg); error(msg);
end

volid = deblank(volid);
len = length(volid);

% Cannot be MINC if less than 5 characters %
if(len < 5) isminc = 0; return; end

last4 = volid(len-3:len);
if(strcmp(last4,'.mnc')) isminc = 1; return; end

if(len < 8) isminc = 0; return; end

last7 = volid(len-6:len);
if(strcmp(last7,'.mnc.gz')) isminc = 1; return; end

isminc = 0;

return;
