% fast_selfreqavg.m - selective frequency averaging
%
% $Id: fast_selfreqavg.m,v 1.3 2003/08/02 01:01:01 greve Exp $
%
% Things to do:
%  1. Save beta, var, and X
%  2. Write ffx combiner
%  3. Save F, ces, cesvar, cespct, cesvarpct, resvar, resstd
%  4. Make compatible with retinotopy
%  5. Wrapper (requires sess)
%  6. TPX
%  7. Whiten/Mask
%  8. SliceTiming Correction
%  9. Run sign reversal
%  10. Global Delay
%  11. Interface with raw twf plot 
%  14. Skirt Nuisance
%  15. nmaxextreg
%  16. Orthog Nuisance
%  17. Implement ER: TR, parname, PSDWin(TER), gammafit
%           Fourier: TR, ncycles
%  18. Nyquist regressor
%  19. Synthesize

% Things done
%  12. PerRun
%  13. JKRun - done

tic;

sesslist = '';
%sesslist = strvcat(sesslist,'/space/greve/2/users/greve/birn-pilot/mgh-dng-22');
%sesslist = strvcat(sesslist,'/home/greve/sg1/mgh-dng-22');
sesslist = strvcat(sesslist,'/space/greve/2/users/greve/birn-pilot/duke-test2');

TR = 3; % Will be needed for tpx
fsd = 'bold';
ananame = 'sm-per';
runlistfile = 'sm.rlf';
%ananame = 'bh-per';
%runlistfile = 'bh.rlf';
conname = 'omnibus';
funcstem = 'fmcsm5';
inorm = 1;
inormtarg = 1000;
ncycles = 8;
nharmonics = 2; % plus 1 for fundamental 
polyfit    = 2;
extregstem = '';
extregstem = 'mcextreg';
%runlistfile = '';
phsigthresh = 2;
dojkrun = 0;
doperrun = 1;
condXthresh = 10e5;

nSess = size(sesslist,1);
nTask = 2*(nharmonics+1);

if(dojkrun & doperrun)
  fprintf('ERROR: cannot specify jkrun and perrun\n');
  return;
end
perrun = [];
jkrun = [];

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
    
  % Loop over external runs (if needed) %
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
    fprintf('Creating design matrix (%g)\n',toc);

    % Load all the information for the design matrix %
    X = [];
    nNuis = 0;
    nFramesTot = 0;
    MeanVal = zeros(nruns,1);
    for nthrun = 1:nruns

      % Get the number of frames for the nth run %
      funcpath = sprintf('%s/%s/%s/%s',sess,fsd,...
			 runlist(nthrun,:),funcstem);
      [nrows ncols nframes fs nslices e bext] = ...
	  fmri_bfiledim(funcpath);
      if(isempty(nrows))
	fprintf('ERROR: reading volume %s\n',funcpath);
	return;
      end
      DM(nthrun).nframes = nframes;
      
      % Task-related component - could do ER or Fourier
      Xtask = [];
      fundamental = 1/(nframes/ncycles);
      t = [0:nframes-1]';
      for nthharmonic = 0:nharmonics
	freq = fundamental*(nthharmonic+1);
	xc = cos(2*pi*t*freq);
	xs = sin(2*pi*t*freq);
	Xtask = [Xtask xc xs];
      end
      DM(nthrun).task = Xtask;
      
      % Nuisance Regressors - poly drift and external reg%
      % Could have multiple extreg here %
      Xpoly = fast_polytrendmtx(1,nframes,1,polyfit);
      extreg = [];
      if(~isempty(extregstem))
	extregpath = sprintf('%s/%s/%s/%s',sess,fsd,...
			     runlist(nthrun,:),extregstem);
	extreg = fmri_ldbvolume(extregpath);
	if(isempty(extreg))
	  fprintf('ERROR: could not load %s\n',extregstem);
	  return;
	end
	if(size(extreg,3)~=1) extreg = squeeze(extreg)'; %'
	else                  extreg = squeeze(extreg);
	end
	% Normalize External Regressor %
	extreg = extreg - repmat(mean(extreg), [nframes 1]);
	extreg = extreg./repmat(std(extreg), [nframes 1]);
      end
      Xnuis = [Xpoly extreg];
      % Orthog Nuisance wrt Self Here - need to keep mean %
      % Orthog Nuisance wrt Task Here - need to keep mean some how %
      DM(nthrun).nuis = Xnuis;
      nNuis = nNuis + size(Xnuis,2);
      nFramesTot = nFramesTot + nframes;
    end % collecting data across runs
    
    % Now create the full design matrix %
    X = [];
    nNuisSum = 0;
    for nthrun = 1:nruns
      nf = DM(nthrun).nframes;
      nNuisRun = size(DM(nthrun).nuis,2);
      npost = nNuis - nNuisSum - nNuisRun;
      z1 = zeros(nf,nNuisSum);
      z2 = zeros(nf,npost);
      Xrun = [DM(nthrun).task z1 DM(nthrun).nuis z2];
      X = [X; Xrun];
      nNuisSum = nNuisSum + size(DM(nthrun).nuis,2);
    end
    %----------- Done Creating design matrix ----------------------%

    XtX = X'*X;
    condX = cond(XtX);
    fprintf('Design Condition: %g\n',condX);
    if(condX > condXthresh)
      fprintf('ERROR: design is badly conditioned\n');
      return;
    end
    d = diag(inv(XtX));
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
	meanvalfile = sprintf('%s.meanval',funcpath);
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
    for nthslice = 0:nslices-1
      %fprintf('nth slice = %d (%g)\n',nthslice,toc);
      %fprintf('%2d (%g) ',nthslice,toc);
      fprintf('%2d ',nthslice);
      if(rem(nthslice,5)==0) fprintf('\n'); end
      
      % Load data for all runs %
      y = [];
      for nthrun = 1:nruns
	funcpath = sprintf('%s/%s/%s/%s',...
			   sess,fsd,runlist(nthrun,:),funcstem);
	[frun tmp1 tmp2 bhdrstr] = fast_ldbslice(funcpath,nthslice);
	if(isempty(frun))
	  fprintf('ERROR: reading volume %s\n',funcpath);
	  return;
	end
	if(inorm) frun = frun*(inormtarg/MeanVal(nthrun)); end
	nframes = size(frun,3);
	frun = reshape(frun,[nvslice nframes])';
	% Multiply frun by Wrun %
	y = [y; frun];
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

      % Save  (if last whitening loop)%
      hoffpath = sprintf('%s/h-offset',anapath);
      fast_svbslice(reshape(ymn,[nrows ncols]),hoffpath,nthslice,'',bhdrstr);
      
      fsigpath = sprintf('%s/fsig',conpath);
      tmp = -log10(Fsig) .* sign(beta(2,:));
      fast_svbslice(reshape(tmp,[nrows ncols]),fsigpath,nthslice,'',bhdrstr);
    
      ph = 180*atan2(beta(2,:),beta(1,:))/pi;
      indz = find(-log10(Fsig) < phsigthresh);
      ph(indz) = 0; % zero out sub-thresh
      phpath = sprintf('%s/phase',conpath);
      fast_svbslice(reshape(ph,[nrows ncols]),phpath,nthslice,'',bhdrstr);
    end % slice

    % If whiten and not last whitening loop
    %  norm and fix acf, compute W for each run, 

    % End whitening loop here 
  
    fprintf('\n\n');
  end % external run loop 

  fprintf('\n\n');
end % session





