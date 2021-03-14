% fmri_tavgslice
% Computes the temporal average at each voxel over multiple runs.
%
%
%


%
% fmri_tavgslice.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
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

