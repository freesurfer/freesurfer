% fast_selfreqavg.m - selective frequency averaging
%
% $Id: fast_selfreqavg.m,v 1.1 2003/07/25 02:28:41 greve Exp $
%
% Things to do:
%  1. Save beta, var, and X
%  2. Write ffx combiner
%  3. Save F, cespct, ces, resvar, resstd
%  4. Make compatible with retinotopy
%  5. Wrapper (requires sess)
%  6. TPX
%  7. Whiten/Mask
%  8. SliceTiming Correction
%  9. Run sign reversal
%  10. Global Delay
%  11. Interface with raw twf plot 
%  12. PerRun/JKRun
%  13. Skirt Nuisance

tic;

sesslist = '';
sesslist = strvcat(sesslist,'/space/greve/2/users/greve/birn-pilot/mgh-dng-22');

TR = 3; % Will be needed for tpx
fsd = 'bold';
ananame = 'sfatst';
conname = 'omnibus';
funcstem = 'fmcsm5';
inorm = 1;
inormtarg = 1000;
ncycles = 8;
nharmonics = 3; % plus 1 for fundamental 
polyfit    = 2;
extregstem = '';
extregstem = 'mcextreg';
%runlistfile = 'bh.rlf';
runlistfile = '';
phsigthresh = 2;

nsess = size(sesslist,1);
ntask = 2*(nharmonics+1);

%---------------- All the sessions -------------------------%
for nthsess = 1:nsess

  sess = deblank(sesslist(nthsess,:));  
  fprintf('nthsess = %d, sess = %s (%g)\n',nthsess,sess,toc);
  
  fsdpath = sprintf('%s/%s',sess,fsd);
  runlist = fast_runlist(fsdpath,runlistfile);
  if(isempty(runlist))
    fprintf('ERROR: could not get run list from %s\n',fsdpath);
    return;
  end
  nruns = size(runlist,1);

  % Create the analysis and contrast directories %
  anapath = sprintf('%s/%s',fsdpath,ananame);
  mkdirp(anapath);
  conpath = sprintf('%s/%s/%s',fsdpath,ananame,conname);
  mkdirp(conpath);
  
  % This is needed to get the number of slices %
  funcpath0 = sprintf('%s/%s/%s/%s',sess,fsd,runlist(1,:),funcstem);
  [nrows ncols nframes fs nslices endian bext] = fmri_bfiledim(funcpath0);
  if(isempty(nrows))
    fprintf('ERROR: reading volume %s\n',funcpath0);
    return;
  end
  nvslice = nrows*ncols;

  % --------- Process each slice separately ---------- %
  for nthslice = 0:nslices-1
    fprintf('nthsslice = %d (%g)\n',nthslice,toc);
    
    % Load all the information for the design matrix %
    X = [];
    nNuis = 0;
    nFramesTot = 0;
    MeanVal = zeros(nruns,1);
    for nthrun = 1:nruns
      
      % Get the number of frames for the nth run %
      funcpath = sprintf('%s/%s/%s/%s',sess,fsd,...
			 runlist(nthrun,:),funcstem);
      [nrows ncols nframes fs nslices e bext] = fmri_bfiledim(funcpath);
      if(isempty(nrows))
	fprintf('ERROR: reading volume %s\n',funcpath);
	return;
      end
      DM(nthrun).nframes = nframes;

      % Task-related component
      Xtask = [];
      fundamental = 1/(nframes/ncycles);
      t = [0:nframes-1]';
      for nthharmonic = 0:nharmonics
	freq = fundamental*(nthharmonic+1);
	xc = cos(2*pi*t*freq);
	xs = sin(2*pi*t*freq);
	%xc = xc - mean(xc);
	%xc = xc/max(xc);
	%xs = xs - mean(xs);
	%xs = xs/max(xs);
	Xtask = [Xtask xc xs];
      end
      DM(nthrun).task = Xtask;

      % Nuisance Regressors - poly drift and external reg%
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
      DM(nthrun).nuis = Xnuis;
      nNuis = nNuis + size(Xnuis,2);
      nFramesTot = nFramesTot + nframes;
    
      % Read the MeanVal here %
      if(inorm)
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
      
    end % run

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
    
    C = [eye(ntask) zeros(ntask,size(X,2)-ntask)];

    % Load the data %
    y = [];
    for nthrun = 1:nruns
      
      funcpath = sprintf('%s/%s/%s/%s',sess,fsd,runlist(nthrun,:),funcstem);
      [frun tmp1 tmp2 bhdrstr] = fast_ldbslice(funcpath,nthslice);
      if(isempty(frun))
	fprintf('ERROR: reading volume %s\n',funcpath);
	return;
      end
      if(inorm) frun = frun*(inormtarg/MeanVal(nthrun)); end
      nframes = size(frun,3);
      frun = reshape(frun,[nvslice nframes])';
      y = [y; frun];
    
    end % run

    % Analyze %
    ymn = mean(y);
    [beta rvar vdof] = fast_glmfit(y,X);
    [F, Fsig, ces] = fast_fratio(beta,X,rvar,C);

    % Save %
    hoffpath = sprintf('%s/h-offset',anapath);
    fast_svbslice(reshape(ymn,[nrows ncols]),hoffpath,nthslice,'',bhdrstr);

    fsigpath = sprintf('%s/fsig',conpath);
    tmp = -log10(Fsig) .* sign(beta(2,:));
    fast_svbslice(reshape(tmp,[nrows ncols]),fsigpath,nthslice,'',bhdrstr);
    
    ph = 180*atan2(beta(2,:),beta(1,:))/pi;
    indz = find(-log10(Fsig) < phsigthresh);
    ph(indz) = 0;
    phpath = sprintf('%s/phase',conpath);
    fast_svbslice(reshape(ph,[nrows ncols]),phpath,nthslice,'',bhdrstr);
    
    
    %return
    
  end % slice
  
  
end % session





