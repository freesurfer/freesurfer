% fast_swfflac_nbhd_sess.m
tic;

% The following variables must be defined:
% flacfile = '~/links/sg1/xval/flac/sem_assoc.flac';
% sesspath = '~/links/sg1/xval/dng';
% outfspec  = 'fmcsm5-swf-bh';
% contrast
% gscaleuse = 1;
% alpha = 0.5;
% sncor = 0; % Assume(1) or don't (0) that the noise and signal are
% spatially correlated
% pthresh = -log10(sigthreshold) regularization
% ytikreg = .1 % Tikhonov regularization parameter
% synthtarg = 0; % Synth target data with WGN
% okfile = '/tmp/ok.ok';
%
% These may or may not be defined
% rcslist = [49 34 14 24 50 14]; % 1-based
% targflacfile = '~/links/sg1/xval/flac/rest.flac';
% 

fast_swfflac_nbhd_sess_ver = '$Id';

SynthSeed = round(sum(100*clock)); 
randn('state',SynthSeed); 
rand('state',SynthSeed); 

fprintf('%s\n',fast_swfflac_nbhd_sess_ver);
fprintf('contrast = %s\n',contrast);
fprintf('alpha = %g\n',alpha);
fprintf('ytkreg = %g\n',ytikreg);
fprintf('sncor = %d\n',sncor);
fprintf('gscaleuse = %d\n',gscaleuse);
fprintf('synthtarg = %g\n',synthtarg);
fprintf('SynthSeed = %g\n',SynthSeed);


% Delete the okfile, if it exists
tmp = sprintf('rm -f %s',okfile);
unix(tmp);

% Handle the row,col,slice list. 1-based
if(~exist('rcslist','var')) rcslist = []; end
if(~isempty(rcslist))
  nlist = length(rcslist)/3;
  rcslist = reshape(rcslist',[3 nlist])';
end

% Load the flac
flac = fast_ldflac(flacfile);
if(isempty(flac)) return; end

% Set the session and do a customize to get the run list
flac.sess = sesspath;
flac.nthrun  = 1;
flac = flac_customize(flac);
if(isempty(flac)) return; end
nruns = size(flac.runlist,1);

fprintf('Starting session %s (%6.1f)\n',flac.sess,toc);    
fprintf('flac = %s, nruns = %d (%g)\n',flac.name,nruns,toc);    
fprintf('outfspec = %s\n',outfspec);

contrastind = flac_conindex(contrast,flac);
if(isempty(contrastind))
  fprintf('ERROR: contrast %s not found in flac\n',contrast);
  return;
end

% Handle when the target flac is diff than source
if(~exist('targflacfile','var')) targflacfile = []; end
if(~isempty(targflacfile))
  targflac = fast_ldflac(targflacfile);
  if(isempty(targflac)) return; end
  targflac.sess = sesspath;
  targflac.nthrun  = 1;
  targflac = flac_customize(targflac);
  if(isempty(targflac)) return; end
  if(size(targflac.runlist,1) ~= nruns)
    fprintf('ERROR: target flac has diff no of runs\n');
    return;
  end
  fprintf('Target flac is %s\n',targflac.name);
else
  targflac = flac;
end

% Load the mask
mstem = sprintf('%s/%s/masks/%s',flac.sess,flac.fsd,flac.mask);
mask = MRIread(mstem);
if(isempty(mask)) return; end
indmask = find(mask.vol);
Nvmask = length(indmask);
indnotmask = find(~mask.vol);
nv = prod(mask.volsize);
volmaskind = zeros(size(mask.vol));
volmaskind(indmask) = 1:length(indmask);

% Set up rcs_delta
rcs_delta = zeros(27,3);
nth = 1;
for r = -1:1
  for c = -1:1
    for s = -1:1
      rcs_delta(nth,:) = [r c s];
      nth = nth + 1;
    end
  end
end




% Go through each run
for jthrun = 1:nruns
  fprintf('\n\n');
  fprintf('Processing jth run %d (%6.1f)\n',jthrun,toc);
  jflac = targflac;
  jflac.nthrun = jthrun;
  jflac = flac_customize(jflac);
  if(isempty(jflac)) return; end
  jflac = flac_desmat(jflac);
  if(isempty(jflac)) return; end

  indtask = flac_taskregind(jflac);
  jX = jflac.X;

  % Go through each run but leave out the jth
  rjk = [];
  sjk = [];
  for kthrun = 1:nruns
    if(kthrun == jthrun) continue; end
    fprintf('  Processing kth run %d (%6.1f)\n',kthrun,toc);
    kflac = flac;
    kflac.nthrun = kthrun;
    kflac = flac_customize(kflac);
    if(isempty(kflac)) return; end

    % Load the residuals
    fprintf('      Loading residuals (%6.1f)\n',toc);
    rstem = sprintf('%s/%s/fla/%s/%s/res',flac.sess,flac.fsd,flac.name,...
		    flac.runlist(kflac.nthrun,:));
    rrun = MRIread(rstem);
    if(isempty(rrun)) return; end
    rrun.vol = fast_vol2mat(rrun.vol);
    rrun.vol = rrun.vol(:,indmask);
    rjk = [rjk; rrun.vol];
    clear rrun;

    % Load the betas
    fprintf('      Loading betas (%6.1f)\n',toc);
    bstem = sprintf('%s/%s/fla/%s/%s/beta',flac.sess,flac.fsd,flac.name,...
		    flac.runlist(kflac.nthrun,:));
    beta = MRIread(bstem);
    if(isempty(beta)) return; end
    beta.vol = fast_vol2mat(beta.vol);
    beta.vol = beta.vol(:,indmask);

    % Project into the contrast subspace
    C = flac.con(contrastind).C;
    srun = (jX*C'*C)*beta.vol;

    sjk = [sjk; srun];
    clear srun beta;
    
  end % Loop over kthrun
  
  % Make sure the dimensions are consistent
  ntpmin = min(size(sjk,1),size(rjk,1));
  sjk = sjk(1:ntpmin,:);
  rjk = rjk(1:ntpmin,:);
  
  % Compute the expected observable
  yjk = sjk + rjk;
  clear rjk;  
  nframesjk = size(yjk,1);
  
  % Load target data
  fprintf('  Loading target, jth run = %d (%6.1f)\n',jthrun,toc);
  fstem = sprintf('%s/%s/%s/%s',jflac.sess,jflac.fsd,...
		  jflac.runlist(jflac.nthrun,:),jflac.funcstem);
  fprintf('  fstem = %s\n',fstem);
  if(~synthtarg)
    y = MRIread(fstem);
    if(isempty(y)) return; end
  else
    fprintf('  Sythesizing target with WGN\n');
    y = MRIread(fstem,1);
    if(isempty(y)) return; end
    y.vol = randn(y.height,y.width,y.depth,y.nframes);
  end
  y.vol = fast_vol2mat(y.vol);
  Nv = size(y.vol,2);
  ysize = [y.height y.width y.depth];
  
  % Loop through each voxel in the mask
  yswf = zeros(y.nframes,Nvmask);
  for nthind = 1:Nvmask
    if(rem(nthind,10000)==0)
      fprintf('%d/%d %g  (%g)\n',nthind,Nvmask,nthind/Nvmask,toc);
    end
    
    % Create a neighborhood around the target voxel
    indv0 = indmask(nthind);
    [r0 c0 s0] = ind2sub(ysize,indv0);
    r = r0 + rcs_delta(:,1);
    c = c0 + rcs_delta(:,2);
    s = s0 + rcs_delta(:,3);

    % Get the indices into the mask
    indm = sub2indmask(r,c,s,volmaskind,indmask);

    % Extract the time courses 
    yjknbhd = yjk(:,indm);
    sjknbhd = sjk(:,nthind);
    G = inv(yjknbhd'*yjknbhd)*yjknbhd'*sjknbhd;
    yraw = y.vol(:,indmask(indm));
    yswf(:,nthind) = yraw*G;
    
  end
  
  % Set mean to be the same as the original
  ymn = mean(y.vol(:,indmask));  
  yswfmn = mean(yswf);
  yswf = yswf + repmat((ymn-yswfmn),[y.nframes 1]);

  % Prepare to save by demasking and reshaping
  yhat = y;
  yhat.vol = zeros(y.nframes,y.height*y.width*y.depth);
  yhat.vol(:,indnotmask) = y.vol(:,indnotmask);
  yhat.vol(:,indmask) = yswf;
  yhat.vol = fast_mat2vol(yhat.vol, [y.height y.width y.depth]);
  
  fprintf('  Saving (%6.1f)\n',toc);    
  yhatfspec = sprintf('%s/%s/%s/%s%s',jflac.sess,jflac.fsd,...
		      jflac.runlist(jflac.nthrun,:),outfspec,...
		      jflac.formatext);
  fprintf('  Saving to %s (%6.1f)\n',yhatfspec,toc);    
  MRIwrite(yhat,yhatfspec);

  logfile = sprintf('%s/%s/%s/%s.log',jflac.sess,jflac.fsd,...
		      jflac.runlist(jflac.nthrun,:),outfspec);
  fp = fopen(logfile,'w');
  fprintf(fp,'%s\n',fast_swfflac_nbhd_sess_ver);
  fprintf(fp,'contrast  = %s\n',contrast);
  fprintf(fp,'alpha     = %g\n',alpha);
  fprintf(fp,'ytkreg    = %g\n',ytikreg);
  fprintf(fp,'sncor     = %d\n',sncor);
  fprintf(fp,'synthtarg = %g\n',synthtarg);
  fprintf(fp,'SynthSeed = %g\n',SynthSeed);
  fclose(fp);

  clear y yswf yhat;

end
fprintf('\n\n');
  
fmri_touch(okfile);
fprintf('Done for session %s (%6.1f)\n',flac.sess,toc);    










