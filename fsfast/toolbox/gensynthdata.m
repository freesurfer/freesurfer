

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