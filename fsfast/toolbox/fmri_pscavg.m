% fmri_pscavg.m
% Purpose: combines percent signal changes of HDRs across sessions.
%   Written as a one-shot program for Dr. Capalan's visit.
% Variables: InputFiles, OutputFile
%
% $Id: fmri_pscavg.m,v 1.1 2003/03/04 20:47:40 greve Exp $



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
