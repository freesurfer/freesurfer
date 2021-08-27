% fast_flafit0.m
%
%
%
% 1. Implement fourier/Run sign reversal/Skirt Nuisance
% 2. Save in appropriate fourier format/Interface with raw twf plot 
% 3. tpx
%  7. Whiten/Mask
%  8. SliceTiming Correction
%  9. nskip
%  10. Global Delay
%  16. Orthog Nuisance
%  18. Nyquist regressor
%  19. Synthesize


%
% fast_flafit0.m
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

%  15. extreg/nmaxextreg
%  20. Multiple conditions
% 2. Save in appropriate sxa format/Interface with raw twf plot 

tic;

sesslist = '';
sesslist = strvcat(sesslist,'/home/greve/sg1/bert-functional');

TR = 2; % Will be needed for tpx
fsd = 'bold';

ananame = 'flatst4-fir';
runlistfile = '';
conname = 'omnibus';
funcstem = 'f';
inorm = 1;
RescaleTarget = 1000;

TimeWindow = 28;
tPreStim =  4;
TER = TR;

Nc = 1;
%parname = 'sem_assoc.par';
parname = 'par01.dat';

hrfmodel = 'fir';
gfDelta = 2.25;
gfTau   = 1.25;

nskip = 0;
tpxlist = [];

%ncycles = 8;
nharmonics = 2; % plus 1 for fundamental 
phsigthresh = 2;

nyqreg = 0;
pfOrder    = 2;
nExtReg    = 3;
extregstem = '';
%extregstem = 'mcextreg';
%runlistfile = '';

dojkrun = 0;
doperrun = 0;
condXthresh = 10e5;

nSess = size(sesslist,1);

if(dojkrun & doperrun)
  fprintf('ERROR: cannot specify jkrun and perrun\n');
  return;
end
mthruntype = 0;
if(doperrun) mthruntype = 1; end
if(dojkrun)  mthruntype = 2; end
perrun = [];
jkrun = [];


%------------- Set up the FLA ------------------------------------%
psdmin = -tPreStim;
dpsd = TER;
psdmax = TimeWindow-tPreStim;
flacfg = fast_flacfg_struct;
flacfg.flaname = ananame;
flacfg.TR = TR;
if(strcmp(hrfmodel,'fir') | strcmp(hrfmodel,'gamma'))
  for c = 1:Nc
    switch(hrfmodel)
     case {'fir'}
      fxline = sprintf('effect fixed cond%d fir %d %g %g %g 0',...
		       c,c,psdmin,dpsd,psdmax);
      Navgs_per_cond = fast_psdwin([psdmin dpsd psdmax],'nirfpsdwin');
      nTask = Nc*Navgs_per_cond;
      gfDelta = [];
      gfTau   = [];
     case {'gamma'}
      fxline = sprintf('effect fixed cond%d gamma %d %g %g %g 0 0 %g %g 0',...
		       c,c,psdmin,dpsd,psdmax,gfDelta,gfTau);
      Navgs_per_cond = 1;
      nTask = Nc;
    end
    flacfg.fxlist(c).fx = fast_fxcfg('parseline',fxline);
  end
% else fourier
end
flacfg.fxlist(Nc+1).fx = fast_fxcfg('parseline',...
   [sprintf('effect random drift polynomial %d',pfOrder)]);
if(~isempty(extregstem))
  flacfg.fxlist(Nc+2).fx = fast_fxcfg('parseline',...
   sprintf('effect random mcextreg extreg %s %d',extregstem, nExtReg));
end

flacfg.fsd = fsd;
flacfg.funcstem = funcstem;
flacfg.evschfname = parname;

flacfg.sesscfg = fast_sesscfg_struct;
ok = fast_fxcfg('checkermid',flacfg);
%-------------------------------------------------------------------%

%---------------- All the sessions -------------------------%
for nthsess = 1:nSess

  sess = deblank(sesslist(nthsess,:));  
  fprintf('nthsess = %d, sess = %s (%g)\n',nthsess,sess,toc);
  
  fsdpath = sprintf('%s/%s',sess,fsd);
  runlist0 = fast_runlist(fsdpath,runlistfile);
  if(isempty(runlist0))
    fprintf('ERROR: could not get run list from %s\n',fsdpath);
    return;
  end
  nruns0 = size(runlist0,1);

  % Determine max for external run loop %
  if(dojkrun | doperrun) nextruns = nruns0;
  else                   nextruns = 1;
  end
    
  % Loop over external runs (if needed) ---------------- %
  for nthextrun = 1:nextruns

    % Create the analysis and omnibus contrast directories %
    if(~doperrun & ~dojkrun)
      anapath = sprintf('%s/%s',fsdpath,ananame);
      runlist = runlist0;
    elseif(doperrun)
      perrun = nthextrun;
      runid = runlist0(perrun,:);
      anapath = sprintf('%s/%s-%s',fsdpath,ananame,runid);
      fprintf('Per Run Loop: nthextrun = %d (%g)\n',nthextrun,toc);
      runlist = runlist0(nthextrun,:);
    elseif(dojkrun)
      jkrun = nthextrun;
      runid = runlist0(jkrun,:);
      anapath = sprintf('%s/%s-jk%s',fsdpath,ananame,runid);
      fprintf('JK Run Loop: nthextrun = %d (%g)\n',nthextrun,toc);
      ind = find([1:nruns0] ~= jkrun);
      runlist = runlist0(ind,:);
    end
    conpath = sprintf('%s/%s',anapath,conname);
    mkdirp(anapath);
    mkdirp(conpath);
    nruns = size(runlist,1);
    
    fprintf('Run List: ');
    for nthrun = 1:nruns
      fprintf('%s ',runlist(nthrun,:));
    end
    fprintf('\n');
    
    % This is needed to get the number of slices %
    funcpath0 = sprintf('%s/%s/%s/%s',sess,fsd,runlist(1,:),funcstem);
    [nrows ncols nframes fs nslices endian bext] = ...
	fmri_bfiledim(funcpath0);
    if(isempty(nrows))
      fprintf('ERROR: reading volume %s\n',funcpath0);
      return;
    end
    nvslice = nrows*ncols;
    
    %----------- Create design matrix ----------------------%
    % Move inside slice loop to do slice delay %
    % set flacfg.tDelay to be slice timing delay %
    fprintf('Creating design matrix (%g)\n',toc);
    flacfg.sesspath = sess;
    
    ok = fast_fxcfg('checkermid',flacfg);
    flacfg = fast_fxcfg('loadsesscfg',flacfg);
    if(isempty(flacfg)) return; end
    X = fast_fla_desmat(flacfg,nthextrun,mthruntype);
    if(isempty(X)) return; end
    %-------------------------------------------------------%
    
    XtX = X'*X;
    condX = cond(XtX);
    fprintf('Design Condition: %g\n',condX);
    if(condX > condXthresh)
      fprintf('ERROR: design is badly conditioned\n');
      return;
    end
    iXtX = inv(XtX);
    hCovMtx = iXtX(1:nTask,1:nTask);
    d = diag(iXtX);
    d = d(1:nTask);
    eff = 1/sum(d);
    fprintf('Design Efficiency: %g\n',eff);
    vrf = 1./d;
    vrfmn = mean(vrf);
    fprintf('VRF: Mean = %g, Min = %g, Max = %g\n',...
	    mean(vrf),min(vrf),max(vrf));

    % Total number of regressors 
    nBeta = size(X,2); 
    
    % Omnibus Contrast %
    C = [eye(nTask) zeros(nTask,nBeta-nTask)]; 
    
    % -------- Read the MeanVal (make sure all exist) ------------ %
    if(inorm)
      for nthrun = 1:nruns
	runid = runlist(nthrun,:);
	meanvalfile = sprintf('%s/%s/%s/%s.meanval',sess,fsd,runid,funcstem);
	fid = fopen(meanvalfile,'r');
	if(fid == -1)
	  fprintf('ERROR: cannot open %s\n',meanvalfile);
	  return;
	end
	MeanVal(nthrun) = fscanf(fid,'%f',1);
	fclose(fid);
	%fprintf('MeanVal = %g\n',MeanVal(nthrun));
      end
    end

    fprintf('\n\n');
    
    % Start whitening loop here (can be more than 2)
    % --------- Process each slice separately ---------- %
    fprintf('Processing data (%g)\n',toc);
    fprintf('Slice: ');
    for nthslice = 0:nslices-1
      fprintf('%2d ',nthslice);
      if(rem(nthslice,10)==9) fprintf('\n       '); end
      
      % Load data for all runs %
      y = [];
      for nthrun = 1:nruns
	funcpath = sprintf('%s/%s/%s/%s',...
			   sess,fsd,runlist(nthrun,:),funcstem);
	[yrun mristruct] = fast_ldbslice(funcpath,nthslice);
	if(isempty(yrun))
	  fprintf('ERROR: reading volume %s\n',funcpath);
	  return;
	end
	if(inorm) 
	  RescaleFactor = (RescaleTarget/MeanVal(nthrun));
	  yrun = yrun*RescaleFactor; 
	else
	  RescaleFactor = 1;
	end
	nframes = size(yrun,3);
	yrun = reshape(yrun,[nvslice nframes])';
	% Multiply yrun by Wrun %
	y = [y; yrun];
      end 
      
      ymn = mean(y);

      % Analyze %
      % Multiply X by Wall %
      [beta rvar vdof] = fast_glmfit(y,X);
      % If whiten and not last whitening loop
      %  get residuals
      %  compute acf (unless last whitening loop) 
      %  sum across mask
      
      % Compute Omnibus Contrast (if last whitening loop)%
      % Multiply X by Wall %
      [F, Fsig, ces] = fast_fratio(beta,X,rvar,C);

      if(strcmp(hrfmodel,'fir') | strcmp(hrfmodel,'gamma'))
	% Save  (if last whitening loop)%
	% Offset from 1st Run
	hoffstem = sprintf('%s/h-offset',anapath);
	fast_svbslice(reshape(beta(nTask+1,:),[nrows ncols]),hoffstem,nthslice,'',mristruct);

	hstem = sprintf('%s/h',anapath);
	hsxa = fast_beta2sxa(beta,rvar,Nc,Navgs_per_cond,X);
	hsxa = reshape(hsxa,[Navgs_per_cond*2*(Nc+1) nrows ncols ]);
	hsxa = permute(hsxa,[2 3 1]);
	fast_svbslice(hsxa,hstem,nthslice,'',mristruct);
	
	fsigpath = sprintf('%s/fsig',conpath);
	tmp = -log10(Fsig) .* sign(sum(ces,1));
	fast_svbslice(reshape(tmp,[nrows ncols]),fsigpath,nthslice,'',mristruct);
	
	fpath = sprintf('%s/f',conpath);
	fast_svbslice(reshape(F,[nrows ncols]),fpath,nthslice,'',mristruct);
    
      end
      
      %ph = 180*atan2(beta(2,:),beta(1,:))/pi;
      %indz = find(-log10(Fsig) < phsigthresh);
      %ph(indz) = 0; % zero out sub-thresh
      %phpath = sprintf('%s/phase',conpath);
      %fast_svbslice(reshape(ph,[nrows ncols]),phpath,nthslice,'',mristruct);
    
    end % slice

    % If whiten and not last whitening loop
    %  norm and fix acf, compute W for each run, 

    % End whitening loop here 

    if(strcmp(hrfmodel,'fir') | strcmp(hrfmodel,'gamma'))

      xfile = sprintf('%s/X.mat',anapath);
      Nnnc = Nc;
      Xfinal = X;
      save(xfile,'Xfinal','Nnnc','pfOrder','nExtReg',...
	   'nruns','Navgs_per_cond','TimeWindow','tPreStim','TR','TER',...
	   'gfDelta','gfTau','tpxlist','RescaleFactor','RescaleTarget',...
	   'nyqreg','-v4');
      
      hdname = sprintf('%s/h.dat',anapath);
      hd = fmri_hdrdatstruct;
      hd.TR  = TR;
      hd.TER = TER;
      hd.TimeWindow = TimeWindow;
      hd.TPreStim = tPreStim;
      hd.Nc = Nc+1;
      hd.Nh = Navgs_per_cond;
      hd.Nnnc = Nc;
      hd.DOF= vdof;
      hd.Npercond = zeros(Nc+1,1);;
      hd.Nruns = nruns;
      hd.Ntp = size(X,1);
      hd.Nrows = nrows;
      hd.Ncols = ncols;
      hd.Nskip = nskip;
      hd.DTOrder = pfOrder + 1;
      hd.RescaleFactor = RescaleFactor;
      hd.HanningRadius = 0.0;
      hd.BrainAirSeg = 0;
      hd.GammaFit = strcmp(hrfmodel,'gamma');
      hd.gfDelta  = gfDelta;
      hd.gfTau         = gfTau;
      hd.NullCondId    = 0;
      hd.SumXtX        = XtX(1:nTask,1:nTask);
      hd.nNoiseAC      = 0;
      hd.CondIdMap     = [0:Nc];
      hd.hCovMtx       = hCovMtx;
      hd.WhitenFlag    = 0;
      hd.runlist       = [];
      for n = 1:nruns, hd.runlist(n) = sscanf(runlist(n,:),'%d'); end
      hd.funcstem      = funcstem;
      hd.parname       = parname;
      if(~isempty(extregstem))
	hd.extregstem  = extregstem;
	hd.nextreg  = nExtReg;
	hd.extortho = 0;
      end
      fmri_svdat3(hdname,hd);
    end

    fprintf('\n\n');
  end % external run loop 

  fprintf('\n\n');
end % session

fprintf('done %g\n\n',toc);




