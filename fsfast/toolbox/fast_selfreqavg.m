% fast_selfreqavg.m - selective frequency averaging
%
%
%
% Things to do:   (nuisance)
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


%
% fast_selfreqavg.m
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


if(0)
%d = '/autofs/space/greve_002/users/greve/birn-pilot';
%d = '/autofs/space/greve_002/users/greve/fb-104.2';
topdir = '/autofs/space/greve_002/users/greve/fbirn-hp-fsfast';

sitelist = '';
sitelist = strvcat(sitelist,'min');
sitelist = strvcat(sitelist,'mgh');
sitelist = strvcat(sitelist,'dunc');
sitelist = strvcat(sitelist,'ucsd');
sitelist = strvcat(sitelist,'nm');
sitelist = strvcat(sitelist,'uci');

sesslist = '';
for site = 1:size(sitelist,1)
  siteid = deblank(sitelist(site,:));
  for subj = [1 3:6]
    for visit = 1:2
      % uci does not have bh and rest for 103.2 %
      if(strcmp(siteid,'uci') & subj==3 & visit==2) continue; end
      if(strcmp(siteid,'mgh')==0)
	sessid = sprintf('%s-data/%s-10%d.%d',siteid,siteid,subj,visit);
      else
	sessid = sprintf('mgh-10%d.%d',subj,visit);
      end	
      sesslist  = strvcat(sesslist,sessid);
    end
  end
end

TR = 3; % Will be needed for tpx
fsd = 'bold';
runlistfile = 'rest.rlf';

ananame = 'rest-sm5-per';
funcstem = 'fmcsm5';
inorm = 1;
doperrun = 1;

%ananame = 'sm-mc-per';
%funcstem = 'mcextreg';
%inorm = 0;

%runlistfile = '';
conname = 'omnibus';
inormtarg = 1000;
Tcycle = 30; % Period of cycle in seconds
nharmonics = 1; % plus 1 for fundamental 
polyfit    = 2;
extregstem = '';
%extregstem = 'mcextreg';
phsigthresh = 2;
dojkrun = 0;
condXthresh = 10e5;
svsignal = 0;
sveres    = 0;
end


%---------------------------------------------------%
tic;
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
  fprintf('nthsess = %d/%d, sess = %s (%g)---------------\n',...
	  nthsess,nSess,sess,toc);
  
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

  % Load the mask, if needed (used to create summary) 
  mask = [];
  if(~isempty(maskstem))
    maskpath = sprintf('%s/%s/masks/%s',sess,fsd,maskstem);
    mask = fast_ldbslice(maskpath);
    if(isempty(mask)) return; end
    nmask = length(find(mask));
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
      fprintf(' Per Run Loop: nthextrun = %d (%g)\n',nthextrun,toc);
      fprintf('   Sess = %s\n',sess);
      runlist = runlist0(nthextrun,:);
    elseif(dojkrun)
      jkrun = nthextrun;
      runid = runlist0(jkrun,:);
      anapath = sprintf('%s/%s-jk%s',fsdpath,ananame,runid);
      fprintf('  JK Run Loop: nthextrun = %d (%g)\n',nthextrun,toc);
      ind = find([1:nruns0] ~= jkrun);
      runlist = runlist0(ind,:);
    end
    conpath = sprintf('%s/%s',anapath,conname);
    fundconpath = sprintf('%s/%s',anapath,'fund');
    estsnrpath = sprintf('%s/estsnr',anapath);
    mkdirp(anapath);
    mkdirp(conpath);
    mkdirp(fundconpath);
    mkdirp(estsnrpath);

    if(sveres)
      eresdir = sprintf('%s/eres',anapath);
      mkdirp(eresdir);
    end
    if(svsignal)
      signaldir = sprintf('%s/signal',anapath);
      mkdirp(signaldir);
    end
    
    nruns = size(runlist,1);

    fprintf('   Run List: ');
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
    fprintf('   Creating design matrix (%g)\n',toc);

    % Load all the information for the design matrix %
    X = [];
    nNuis = 0;
    nFramesTot = 0;
    MeanVal = zeros(nruns,1);
    nextregtot = 0;
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
      %fundamental = 1/(nframes/ncycles);
      fundamental = 1/Tcycle;
      t = TR*[0:nframes-1]';
      for nthharmonic = 0:nharmonics
	freq = fundamental*(nthharmonic+1);
	xc = cos(2*pi*t*freq);
	%xc = xc - mean(xc);
	xs = sin(2*pi*t*freq);
	%xs = xs - mean(xs);
	Xtask = [Xtask xc xs];
      end
      DM(nthrun).task = Xtask;
      
      % Nuisance Regressors - poly drift and external reg%
      % Could have multiple extreg here %
      Xpoly = fast_polytrendmtx(1,nframes,1,polyfit);
      DM(nthrun).poly = Xpoly;
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
	if(~isempty(nextreg)) extreg = extreg(:,1:nextreg); end
	nextregtot = nextregtot + size(extreg,2);
	
	% Demean and Normalize External Regressor %
	extreg = extreg - repmat(mean(extreg), [nframes 1]);
	extreg = extreg./repmat(std(extreg), [nframes 1]);
	DM(nthrun).extreg = extreg;
	DM(nthrun).Rextreg = eye(nframes) - extreg*inv(extreg'*extreg)*extreg';
      
	if(extregorthog & 0)
	  fprintf('   Orthognalizing design wrt external regressor\n');
	  Xtask = DM(nthrun).Rextreg*Xtask;
	  DM(nthrun).task = Xtask;
	end
      end

      if(extregorthog)
	Xnuis = Xpoly;
      else
	Xnuis = [Xpoly extreg];
      end

      DM(nthrun).nuis = Xnuis;
      nNuis = nNuis + size(Xnuis,2);
      nFramesTot = nFramesTot + nframes;
    
    end % collecting design info across runs
    
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
    fprintf('   Design Condition: %g\n',condX);
    if(condX > condXthresh)
      fprintf('ERROR: design is badly conditioned\n');
      return;
    end
    d = diag(inv(XtX));
    d = d(1:nTask);
    eff = 1/sum(d);
    fprintf('   Design Efficiency: %g\n',eff);
    vrf = 1./d;
    vrfmn = mean(vrf);
    fprintf('   VRF: Mean = %g, Min = %g, Max = %g\n',...
	    mean(vrf),min(vrf),max(vrf));

    % Total number of regressors 
    nBeta = size(X,2); 
    
    % Omnibus Contrast %
    C = [eye(nTask) zeros(nTask,nBeta-nTask)]; 

    % Fundamental Contrast %
    Cfund = zeros(1,nBeta);
    Cfund(1:2) = 1; % Real and Imag

    % Save info to X.mat
    xmatpath = sprintf('%s/X.mat',anapath);
    save(xmatpath,'X','DM','sess','fsd','funcstem',...
	 'runlist','nTask','C','inorm','inormtarg','Tcycle',...
	 'nharmonics','polyfit','extregstem','nextreg','extregorthog',...
	 'phsigthresh','dojkrun','doperrun','condXthresh');
    
    % -------- Read the MeanVal (make sure all exist) ------------ %
    if(inorm)
      for nthrun = 1:nruns
	runid = runlist(nthrun,:);
	meanvalfile = sprintf('%s/%s/%s/%s.meanval',...
			      sess,fsd,runid,funcstem);
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

    fprintf('\n');
    
    % Start whitening loop here (can be more than 2)
    % --------- Process each slice separately ---------- %
    fprintf('  Processing data (%g)\n   ',toc);
    rvarpctsum = 0;
    ar1sum    = 0;
    nover     = 0;
    cnrsum    = 0;
    magpctsum = 0;
    for nthslice = 0:nslices-1
      %fprintf('nth slice = %d (%g)\n',nthslice,toc);
      %fprintf('%2d (%g) ',nthslice,toc);
      fprintf('%2d ',nthslice);
      if(rem(nthslice+1,15)==0) fprintf('\n   '); end
      
      % Load data for all runs %
      y = [];
      clear nframes_per_run;
      for nthrun = 1:nruns
	funcpath = sprintf('%s/%s/%s/%s',...
			   sess,fsd,runlist(nthrun,:),funcstem);
	[frun mristruct] = fast_ldbslice(funcpath,nthslice);
	if(isempty(frun))
	  fprintf('ERROR: reading volume %s\n',funcpath);
	  return;
	end
	if(inorm) frun = frun*(inormtarg/MeanVal(nthrun)); end
	nframes = size(frun,3);
	nframes_per_run(nthrun) = nframes;
	frun = reshape(frun,[nvslice nframes])';
	if(extregorthog)
	  %fprintf('   Orthognalizing design wrt external regressor\n');
	  frun = DM(nthrun).Rextreg*frun;
	end
	% Multiply frun by Wrun %
	y = [y; frun];
      end 
      
      ymn = mean(y);
      ind = find(ymn == 0);
      ymn(ind) = eps;

      % Analyze %
      % Multiply X by Wall %
      [beta rvar vdof r] = fast_glmfit(y,X);

      % Adjust the DOF if ExtReg were projected out
      if(extregorthog)	vdof = vdof - nextregtot;  end
      
      if(sveres)
	i1 = 1;
	for nthrun = 1:nruns
	  nframes = nframes_per_run(nthrun);
	  i2 = i1 + nframes - 1;
	  erestmp = reshape(r(i1:i2,:)',[nrows ncols nframes]);
	  stem = sprintf('%s/e%03d',eresdir,nthrun);
	  fast_svbslice(erestmp,stem,nthslice,'',mristruct);
	  i1 = i2 + 1;
	end
      end
      if(svsignal)
	i1 = 1;
	for nthrun = 1:nruns
	  nframes = nframes_per_run(nthrun);
	  i2 = i1 + nframes - 1;
	  sigtask = X(i1:i2,1:nTask)*beta(1:nTask,:);
	  sigtask = reshape(sigtask',[nrows ncols nframes]);
	  stem = sprintf('%s/s%03d',signaldir,nthrun);
	  fast_svbslice(sigtask,stem,nthslice,'',mristruct);
	  i1 = i2 + 1;
	end
      end
      
      
      % If whiten and not last whitening loop
      %  get residuals
      %  compute acf (unless last whitening loop) 
      %  sum across mask
      
      % Compute Omnibus Contrast (if last whitening loop)%
      % Multiply X by Wall %
      [F, Fsig, ces] = fast_fratio(beta,X,rvar,C,[],[],vdof);

      % Signifiance of F-test %
      fsigpath = sprintf('%s/fsig',conpath);
      tmp = -log10(Fsig) .* sign(beta(2,:));
      fast_svbslice(reshape(tmp,[nrows ncols]),fsigpath,nthslice,'',mristruct);
    
      % F-test ratio %
      fpath = sprintf('%s/f',conpath);
      fast_svbslice(reshape(F,[nrows ncols]),fpath,nthslice,'',mristruct);
    
      % F Test for the fundamental
      [Ffund, Fsigfund, cesfund] = fast_fratio(beta,X,rvar,Cfund,[],[],vdof);
      fsigpath = sprintf('%s/fsig',fundconpath);
      tmp = -log10(Fsigfund) .* sign(beta(2,:));
      fast_svbslice(reshape(tmp,[nrows ncols]),fsigpath,nthslice,'',mristruct);
      fpath = sprintf('%s/f',fundconpath);
      fast_svbslice(reshape(Ffund,[nrows ncols]),fpath,nthslice,'',mristruct);
    
      % Save results to disk (if last whitening loop)%

      % Mean offset %
      hoffpath = sprintf('%s/h-offset',anapath);
      fast_svbslice(reshape(ymn,[nrows ncols]),hoffpath,nthslice,'',mristruct);
      
      % Beta %
      betapath = sprintf('%s/beta',anapath);
      tmp = reshape(beta',[nrows ncols nBeta ]);
      fast_svbslice(tmp,betapath,nthslice,'',mristruct);
      
      % Residual error variance %
      rvarpath = sprintf('%s/rvar',estsnrpath);
      tmp = reshape(rvar,[nrows ncols]);
      fast_svbslice(tmp,rvarpath,nthslice,'',mristruct);
      
      % Residual error ar1 %
      rar0 = mean(r.*r);
      indz = find(rar0==0);
      rar0(indz) = eps;
      rar1 = mean( r(1:end-1,:) .* r(2:end,:) ) ./ rar0;
      rar1path = sprintf('%s/ar1',estsnrpath);
      tmp = reshape(rar1,[nrows ncols]);
      fast_svbslice(tmp,rar1path,nthslice,'',mristruct);
      
      % Magnitude of fundamental
      magpct = 100*sqrt( beta(1,:).^2 + beta(2,:).^2 ) ./ ymn;
      magpath = sprintf('%s/magpct',conpath);
      tmp = reshape(magpct,[nrows ncols]);
      fast_svbslice(tmp,magpath,nthslice,'',mristruct);
      
      % Phase of fundamental %
      ph = 180*atan2(beta(2,:),beta(1,:))/pi;
      indz = find(-log10(Fsig) < phsigthresh);
      ph(indz) = 0; % zero out sub-thresh
      phpath = sprintf('%s/phase',conpath);
      fast_svbslice(reshape(ph,[nrows ncols]),phpath,nthslice,'',mristruct);
    
      % Create summary table %
      if(~isempty(mask))
	indmask_in  = find(mask(:,:,nthslice+1));
	indmask_out = find(~mask(:,:,nthslice+1));

	rvarpct = 100*rvar./ymn;
	rvarpctsum = rvarpctsum + sum(rvarpct(indmask_in));

	ar1sum = ar1sum + sum(rar1(indmask_in));

	pomni = -log10(Fsig);
	pomni(indmask_out) = 0;
	pthresh = 3;
	indover = find(abs(pomni) > pthresh);
	nover = nover + length(indover);

	cnrsum = cnrsum + sum(F(indover));
	magpctsum = magpctsum + sum(magpct(indover));
      end
        
    end % slice

    % If whiten and not last whitening loop
    %  norm and fix acf, compute W for each run, 

    % End whitening loop here 
  
    if(~isempty(mask))
      rstdpctmn = sqrt(rvarpctsum/nmask);
      ar1mn     = ar1sum/nmask;
      cnrmn     = cnrsum/nover;
      magpctmn  = magpctsum/nover;
      pctover   = 100*nover/nmask;
      outfile = sprintf('%s/table.sum',anapath);
      fprintf('Saving summary table to %s\n',outfile);
      fp = fopen(outfile,'w');
      fprintf(fp,'%5d  %5.2f  %5.2f  %5.2f  %4.2f  %4.2f\n',...
	      nmask,pctover,magpctmn,cnrmn,rstdpctmn,ar1mn);
      fclose(fp);
    end
    
    fprintf('\n\n');
  end % external run loop 

  fprintf('\n\n');
end % session

fprintf('fast_selfreqavg done (%g)\n\n',toc);



