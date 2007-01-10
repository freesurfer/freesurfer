function tTPExclude = fmri_ldtpexcl(TPExclFile)
%
% tTPExclude = ldtpexcl(TPExclFile)
%
% Reads the TP Exclude File.
%
%
%


%
% fmri_ldtpexcl.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
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

tTPExclude = [];

TPExclFile = deblank(TPExclFile);

if( strcmp(TPExclFile,'noexcl')) return; end

fid = fopen(TPExclFile);
if(fid == -1)
  msg = sprintf('Could not open %s',TPExclFile);
  qoe(msg);
  error(msg);
end

tTPExclude = fscanf(fid,'%f');

return;
  
