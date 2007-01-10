

%
% gensynthdata.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
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


TopDir = '/space/raid/5/users/greve/synthdata';
SessionDir = strcat(TopDir,'/WSNR-100c');
ParDir     = strcat(SessionDir,'/par');
SNR = 1;

ucmd = sprintf('mkdir -p %s',ParDir);
unix(ucmd);

nRuns = 1;

for run = 1:nRuns,

  fprintf('-------- Generating Run %d ----------\n',run);
  srs = fmri_synthrunstruct;
  srs.SNR = SNR;
  srs.nPerCond   = [32 32 32 32];

  [srs fSlices] = fmri_synthrun(srs);

  ParFile = sprintf('%s/par%02d.dat',ParDir,run);
  fmri_svpar(srs.Par,ParFile);

  RunDir = sprintf('%s/%02d',SessionDir,run);
  ucmd = sprintf('mkdir -p %s',RunDir);
  unix(ucmd);

  for slice = 0:srs.nSlices-1,
   SliceFile = sprintf('%s/f_%03d.bshort',RunDir,slice);
   fmri_svbfile(fSlices(:,:,:,slice+1),SliceFile);
  end

end
