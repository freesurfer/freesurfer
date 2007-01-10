% fmri_tavgslice
% Computes the temporal average at each voxel over multiple runs.
%
%
%


%
% fmri_tavgslice.m
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

fprintf(1,'\n');

nRuns = size(InputFiles,1);

Ntp = 0;
ysum = 0;
for r = 1:nRuns,
  fprintf(1,'Loading %s\n',InputFiles(r,:));
  y = fmri_ldbfile(InputFiles(r,:));
  ysum = ysum + sum(y,3);
  Ntp = Ntp + size(y,3);
end

ysum = ysum/Ntp;

fprintf(1,'Saving %s\n',tAvgFile);
fmri_svbfile(ysum,tAvgFile);

