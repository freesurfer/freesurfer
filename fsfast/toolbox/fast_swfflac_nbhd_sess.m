% fast_swfflac_nbhd_sess.m


%
% fast_swfflac_nbhd_sess.m
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



fast_swfflac_nbhd_sess_ver = 'fast_swfflac_nbhd_sess.m @FS_VERSION@';

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
fprintf('assumeiid = %d\n',assumeiid);
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
C = flac.con(contrastind).C;

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

% Load reference info
if(0)
clear drill;
drill.brain = mask;
drill.refmask = MRIread('~/links/sg1/xval/dng/bold/masks/samc-rsyn-mask.mgh');
drill.indref = find(drill.refmask.vol);
drill.nref = length(drill.indref);
drill.betaref = MRIread('~/links/sg1/xval/dng/bold/samc/ffx/beta.mgh');
drill.betaref.vol = fast_vol2mat(drill.betaref.vol);
drill.gamref = drill.betaref;
drill.gamref.vol = C*drill.betaref.vol;
drill.F0    = MRIread('~/links/sg1/xval/dng/bold/samc-rsyn/ffx/omnibus/f.mgh');
drill.Fsig0 = MRIread('~/links/sg1/xval/dng/bold/samc-rsyn/ffx/omnibus/fsig.mgh');
drill.Fsnc0    = MRIread('~/links/sg1/xval/dng/bold/samc-rsyn-snc-00-gn/ffx/omnibus/f.mgh');
drill.Fsigsnc0 = MRIread('~/links/sg1/xval/dng/bold/samc-rsyn-snc-00-gn/ffx/omnibus/fsig.mgh');
drill.rdgamssesnc0 = MRIread('~/links/sg1/xval/dng/bold/samc-rsyn-snc-00-gn/ffx/omnibus/rdgamsse.mgh');
end

ntp = flac.ntp;

% Go through each run
for jthrun = 1:nruns
  nthref = 1;
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
  [Fsig F betaffx] = flacffx(flac,contrast,1,jthrun);
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
    rstem = sprintf('%s/%s/%s/%s/res',flac.sess,flac.fsd,flac.name,...
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
      bstem = sprintf('%s/%s/%s/%s/beta',flac.sess,flac.fsd,...
		    flac.name,flac.runlist(kflac.nthrun,:));
      beta = MRIread(bstem);
    else
      beta = betaffx;
    end
    if(isempty(beta)) return; end
    beta.vol = fast_vol2mat(beta.vol);
    beta.vol = beta.vol(:,indmask);

    % Project into the contrast subspace
    srun = (jX*C'*C)*beta.vol;

    sjk = [sjk; srun];
    %clear srun beta;
    
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
  nillcond = 0;
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
    indm0 = indm;
    
    pvaltarg = Fp.vol(indv0);    
    pvalm = Fp.vol(indmask(indm));
    if(fdr < 1 & ~sop)
      % Get the p-values
      nover = length(find(pvalm < pthresh));
      if(pvaltarg > pthresh & nover > 0)
	% Target is not sig but others are
	% Keep ones that are less sig than targ (incl targ)
	indkeep = find(pvalm >= pvaltarg); 
	indm = indm(indkeep);
	nadjust = nadjust + 1;
      else 
	indkeep = [1:length(indm0)]';
      end
    end
    if(sop)
      if(pvaltarg > pthresh)
	% Target is not sig, so keep only insig nbhrs
	indkeep = find(pvalm > pthresh); 
      else
	% Target is sig, so keep only sig nbhrs
	indkeep = find(pvalm <= pthresh); 	
      end
      indm = indm(indkeep);
    end
    nnbrs = length(indm)-1;
    
    % Extract the time courses 
    yraw = y.vol(:,indmask(indm));
    yjknbhd = yjk(:,indm);
    sjknbhd = sjk(:,indm);
    rjknbhd = rjk(:,indm);
    sjktarg = sjk(:,nthind);
    rjktarg = rjk(:,nthind);

    if(~assumeiid)
      if(usesnc) 
	yjknbhd = [yjknbhd yjksnc]; 
	sjknbhd = [sjknbhd sjksnc]; 
      end
      
      % Compute and apply the filter
      yjkscm = yjknbhd'*yjknbhd;
      yjkscm_cond = cond(yjkscm);
      if(yjkscm_cond < 10e10)
	if(sncor)
	  G = inv(yjknbhd'*yjknbhd)*yjknbhd'*sjktarg;
	else
	  G = inv(yjknbhd'*yjknbhd)*sjknbhd'*sjktarg;
	end    
      else
	nillcond = nillcond + 1;
	G = zeros(size(yjknbhd,2),1);
	G(1) = 1;
      end

      if(gnorm)  
	sjktarghat = sjknbhd*G;
	tmp = sjktarghat'*sjktarghat;
	if(tmp ~= 0)
	  Gscale = inv(sjktarghat'*sjktarghat)*sjktarghat'*sjktarg;
	else
	  Gscale = 1;
	end
	G = Gscale*G;
	%G = G/sqrt(sum(G.^2));
      end
      
      if(usesnc) 
	%yraw0 = yraw;
	yraw = [yraw yrawsnc];
      end
      
      %yswfv0 = yraw0*G0;
      yswfv = yraw*G;
    else
      % assume iid
      if(nnbrs > 0)
	snottarg = sjknbhd(:,2:end);
	sm = sum(snottarg.*snottarg);
	indtmp = find(sm==0);
	sm(indtmp) = 10e10;
	G0 = (snottarg'*sjktarg)'./sm;
	G = [1 G0]/size(sjknbhd,2);
	G = G(:);
	%shat0 = snottarg.*repmat(G0,[ntpmin 1]);
	%shat = sjknbhd*G;
	nn = 1:ntpmin;
	%plot(nn,sjktarg,nn,shat);
	yswfv = yraw*G;
      else
	yswfv = yraw;
      end

      if(usesnc & nnbrs > 0)
	rtil = rjknbhd*G;
	Gsnc = inv(rjksnc'*rjksnc)*rjksnc'*rtil;
	Gsnc = Gsnc(:);
	yswfv = yswfv - yrawsnc*Gsnc;
	%rtilhat = rjksnc*Gsnc;
	%d = rtil-rtilhat;
      end
      %if(Fp.vol(indv0)<.01) keyboard; end
    end

    yswf(:,nthind) = yswfv;
    
    
    %if(indv0 == 23151) 
    if(0)
    if(drill.refmask.vol(indv0))
      G0 = zeros(length(indm0),1);
      G0 = G;
      yjknbhd0 = yjk(:,indm0);
      sjknbhd0 = sjk(:,indm0);
      rjknbhd0 = rjk(:,indm0);
      yraw0 = y.vol(:,indmask(indm0));
      drill.indv0(nthref) = indv0;
      drill.v(nthref).indm0 = indm0;
      drill.v(nthref).pvalm0 = pvalm;
      drill.v(nthref).sjknbhd0 = sjknbhd0;
      drill.v(nthref).rjknbhd0 = rjknbhd0;
      %drill.v(nthref).yjknbhd0 = yjknbhd0; % equals s+r
      drill.v(nthref).G0 = G0;
      drill.v(nthref).yraw0 = yraw0;
      drill.yswf(:,nthref)  = yswfv;
      drill.Gscale(nthref) = Gscale;
      nthref = nthref + 1;

      if(0)
      fprintf('nthind = %d, %g\n',nthind,Fp.vol(indmask(nthind)) );
      matfile = sprintf('swf.drill.%d.%d.mat',indv0,jthrun);
      if(str2num(version('-release')) < 14)
	save(matfile,'indv0','jthrun','indmask','indm0','indkeep',...
	     'indm','pthresh','pvaltarg','pvalm','yjknbhd0','sjknbhd0',...
	     'rjknbhd0','yraw0','G','Gscale');
      else
	save(matfile,'indv0','jthrun','indmask','indm0','indkeep',...
	     'indm','pthresh','pvaltarg','pvalm','yjknbhd0','sjknbhd0',...
	     'rjknbhd0','yraw0','G','Gscale','-v6');
      end
      end
    end
    end

  end % loop over mask voxels
  fprintf('nadjust  = %d/%d = %g%%\n',nadjust,Nvmask,100*nadjust/Nvmask);
  fprintf('nillcond = %d/%d = %g%%\n',nillcond,Nvmask,100*nillcond/Nvmask);

  if(0)
  fprintf('Saving drill\n');
  drillmat = sprintf('drill-%d.mat',jthrun);
  if(str2num(version('-release')) < 14)
    save(drillmat,'drill','jX','pthresh');
  else
    save(drillmat,'drill','jthrun','jX','pthresh','-v6');
  end  
  drill.v = []; %clears v but keeps other components
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
  fprintf(fp,'fdr       = %g\n',fdr);
  fprintf(fp,'FDRThreshold is %g\n',pthresh);
  fprintf(fp,'nadjust   = %g (%g%%)\n',nadjust,100*nadjust/Nvmask);
  fprintf(fp,'nillcond  = %g (%g%%)\n',nillcond,100*nillcond/Nvmask);
  fprintf(fp,'useffxbeta = %d\n',useffxbeta);
  fprintf(fp,'usesnc = %d\n',usesnc);
  fprintf(fp,'snconly = %d\n',snconly);
  fprintf(fp,'sncdim = %d\n',sncdim);
  fprintf(fp,'sncor = %d\n',sncor);
  fprintf(fp,'sop = %d\n',sop);
  fprintf('assumeiid = %d\n',assumeiid);
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










