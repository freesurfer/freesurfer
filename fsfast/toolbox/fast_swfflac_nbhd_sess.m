% fast_swfflac_nbhd_sess.m
tic;

% fdr = 1;
% usesnc = 0;
% sncdim = 20;


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

fast_swfflac_nbhd_sess_ver = '$Id: fast_swfflac_nbhd_sess.m,v 1.3 2005/01/23 17:09:21 greve Exp $';

SynthSeed = round(sum(100*clock)); 
randn('state',SynthSeed); 
rand('state',SynthSeed); 

fprintf('%s\n',fast_swfflac_nbhd_sess_ver);
fprintf('contrast = %s\n',contrast);
fprintf('alpha = %g\n',alpha);
fprintf('ytkreg = %g\n',ytikreg);
fprintf('sncor = %d\n',sncor);
fprintf('usesnc = %d\n',usesnc);
fprintf('useffxbeta = %d\n',useffxbeta);
fprintf('snconly = %d\n',snconly);
fprintf('gnorm = %d\n',gnorm);
fprintf('sop = %d\n',sop);
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
if(nbhd_inplane) srange = 0;
else             srange = -1:1;
end
clear rcs_delta;
nth = 1;
for r = -1:1
  for c = -1:1
    for s = srange
      if(r==0 & c == 0 & s == 0) continue; end
      rcs_delta(nth,:) = [r c s];
      nth = nth + 1;
    end
  end
end
if(~snconly)
  fprintf('Neighborhood RCS\n');
  fprintf('%2d %2d %2d\n',rcs_delta');
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

  % Perform FFX with jackknifing
  fprintf('  Computing FFX\n');
  [Fsig F betaffx] = flacffx(flac,contrast,0,jthrun);
  Fp = Fsig;
  Fp.vol = 10.^(-Fsig.vol); % convert to p
  
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
    if(~useffxbeta)
      fprintf('      Loading betas (%6.1f)\n',toc);
      bstem = sprintf('%s/%s/fla/%s/%s/beta',flac.sess,flac.fsd,...
		    flac.name,flac.runlist(kflac.nthrun,:));
      beta = MRIread(bstem);
    else
      beta = betaffx;
    end
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
  nframesjk = size(yjk,1);
  
  % Get voxels for SNC
  if(usesnc)
    fprintf('sncdim = %d\n',sncdim);
    indpas = find(Fp.vol(indmask) > 0.75);
    npas = length(indpas);
    indact = find(Fp.vol(indmask) < .01);
    nact = length(indact);
    fprintf('  npas %d, nact = %d\n',npas,nact);
    fprintf('  Computing SVD (t=%g)\n',toc);
    [Ur Sr Vr] = fast_svd(rjk(:,indpas));
    fprintf('  Done computing SVD (t=%g)\n',toc);
    if(0)
      %[indbest cpvsr] = fast_tspvsrank(rjk(:,indact),Ur(:,1:sncdim),sncdim);
      [indbest cpvsr] = fast_tspvsrank(rjk(:,indact),Ur(:,1:2*sncdim),sncdim);
      fprintf('  IndBest: ');
      fprintf('%d ',indbest(1:sncdim));
      fprintf('  \n');
      fprintf('  CPVS: ');
      fprintf('%6.1f ',cpvsr(1:sncdim));
      fprintf('  \n');
    else
      indbest = 1:sncdim;
    end
    Vr = Vr(:,indbest);
    sjksnc = sjk(:,indpas)*Vr; 
    rjksnc = rjk(:,indpas)*Vr; 
    %yjksnc = yjk(:,indpas)*Vr; 
    yjksnc = rjksnc; % sneaky, sneaky!
    sjksnc = zeros(size(sjksnc)); % sneaky, sneaky!
  end
  %clear rjk;  

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
  
  if(usesnc) yrawsnc = y.vol(:,indmask(indpas))*Vr; end

  pthresh = fast_fdrthresh(Fp.vol(indmask),fdr);
  if(fdr < 1)
    fprintf('  Global FDR threshold is %g (%g)\n',pthresh,-log10(pthresh));
    noverthresh = length(find(Fp.vol(indmask) < pthresh));
    fprintf('    %d (%4.1f%%) over thresh\n',noverthresh,100*noverthresh/Nvmask);
  end

  % Loop through each voxel in the mask
  nadjust = 0;
  yswf = zeros(y.nframes,Nvmask);
  for nthind = 1:Nvmask
    if(rem(nthind,10000)==0)
      fprintf('%d/%d %g  (%g)\n',nthind,Nvmask,nthind/Nvmask,toc);
    end
    
    indv0 = indmask(nthind);
    if(~snconly)
      % Create a neighborhood around the target voxel
      [r0 c0 s0] = ind2sub(ysize,indv0);
      r = r0 + rcs_delta(:,1);
      c = c0 + rcs_delta(:,2);
      s = s0 + rcs_delta(:,3);

      % Get the indices into the mask
      indm = sub2indmask(r,c,s,volmaskind,indmask);
      indm = [nthind; indm]; % put target voxel first
    else
      indm = [nthind]; % only target voxel
    end
    
    if(fdr < 1)
      % Get the p-values
      pvaltarg = Fp.vol(indv0);    
      pvalm = Fp.vol(indmask(indm));
      %pthresh = fast_fdrthresh(pvalm,fdr);% do this globally now
      nover = length(find(pvalm < pthresh));
      if(pvaltarg > pthresh & nover > 0)
	% Target is not sig but others are
	% Keep ones that are less sig than targ (incl targ)
	indkeep = find(pvalm >= pvaltarg); 
	indm = indm(indkeep);
	nadjust = nadjust + 1;
      end
    end
    if(sop)
      % Spill-over protection
      pvaltarg = Fp.vol(indv0);    
      pvalm = Fp.vol(indmask(indm));
      indkeep = find(pvalm >= pvaltarg); 
      indm = indm(indkeep);
    end
    
    % Extract the time courses 
    yjknbhd = yjk(:,indm);
    sjknbhd = sjk(:,indm);
    sjktarg = sjk(:,nthind);
    rjktarg = rjk(:,nthind);

    if(usesnc) 
      yjknbhd0 = yjknbhd;
      sjknbhd0 = sjknbhd;
      yjknbhd2 = [yjknbhd0 rjksnc]; 
      yjknbhd = [yjknbhd yjksnc]; 
      sjknbhd = [sjknbhd sjksnc]; 
    end

    % Compute and apply the filter
    if(sncor)
      %G0 = inv(yjknbhd0'*yjknbhd0)*yjknbhd0'*sjktarg;
      %sjktarg0hat = yjknbhd0*G0;
      G = inv(yjknbhd'*yjknbhd)*yjknbhd'*sjktarg;
    else
      G = inv(yjknbhd'*yjknbhd)*sjknbhd'*sjktarg;
    end    

    if(gnorm)  
      sjktarghat = sjknbhd*G;
      Gscale = inv(sjktarghat'*sjktarghat)*sjktarghat'*sjktarg;
      G = Gscale*G;
      %G = G/sqrt(sum(G.^2));
    end

    yraw = y.vol(:,indmask(indm));
    if(usesnc) 
      %yraw0 = yraw;
      yraw = [yraw yrawsnc];
    end

    %yswfv0 = yraw0*G0;
    yswfv = yraw*G;
    yswf(:,nthind) = yswfv;

    if(0 & Fp.vol(indmask(nthind)) < .001) 
      fprintf('nthind = %d, %g\n',nthind,Fp.vol(indmask(nthind)) );
      keyboard; 
    end

  end % loop over mask voxels
  fprintf('nadjust = %d/%d = %g%%\n',nadjust,Nvmask,100*nadjust/Nvmask);
  
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
  fprintf(fp,'fdr       = %g\n',fdr);
  fprintf(fp,'FDRThreshold is %g\n',pthresh);
  fprintf(fp,'nadjust   = %g (%g%%)\n',nadjust,100*nadjust/Nvmask);
  fprintf(fp,'useffxbeta = %d\n',useffxbeta);
  fprintf(fp,'usesnc = %d\n',usesnc);
  fprintf(fp,'snconly = %d\n',snconly);
  fprintf(fp,'sncdim = %d\n',sncdim);
  fprintf(fp,'sncor = %d\n',sncor);
  fprintf(fp,'sop = %d\n',sop);
  fprintf(fp,'gnorm = %d\n',gnorm);
  fprintf(fp,'nbhd_inplane = %d\n',nbhd_inplane);
  fprintf(fp,'nnbhd = %d\n',size(rcs_delta,1)+1);
  fprintf(fp,'synthtarg = %g\n',synthtarg);
  fprintf(fp,'SynthSeed = %g\n',SynthSeed);
  fclose(fp);

  clear y yswf yhat;

end % Loop over jked rund

fprintf('\n\n');
  
fmri_touch(okfile);
fprintf('Done for session %s (%6.1f)\n',flac.sess,toc);    










