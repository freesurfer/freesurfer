% fmri_pscavg.m
% Purpose: combines percent signal changes of HDRs across sessions.
%   Written as a one-shot program for Dr. Capalan's visit.
% Variables: InputFiles, OutputFile
%
%


%
% fmri_pscavg.m
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



nInputs = size(InputFiles,1);

psc_sum = 0;

for n = 1:nInputs
  fprintf(1,'Opening %s\n',deblank(InputFiles(n,:)));
  fid = fopen(deblank(InputFiles(n,:)),'r');
  if(fid == -1)
    msg = sprintf('Could not open/find %s for reading',deblank(InputFiles(n,:)));
    qoe(msg);error(msg);
  end
  q = fscanf(fid,'%f');
  fclose(fid);

  q = reshape(q, [3+nNNC nHEst])';
  tavgavg  = q(1,3);
  havg     = q(:,4: size(q,2));
  havg = reshape1d(havg);

  psc = [psc; havg'/tavgtavg];

end

psc_avg = mean(psc);
psc_std = std(psc);
