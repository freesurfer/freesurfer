% fmri_tavgslice
% Computes the temporal average at each voxel over multiple runs.
%
%
% $Id: fmri_tavgslice.m,v 1.1 2003/03/04 20:47:40 greve Exp $

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

