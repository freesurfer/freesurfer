function r = fast_selxavg(varargin)
% r = fast_selxavg(varargin)


%
% fast_selxavg.m
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

version = 'fast_selxavg.m @FS_VERSION@';
fprintf(1,'%s\n',version);
r = 1;

%% Print useage if there are no arguments %%
if(nargin == 0)
  print_usage;
  return;
end

%% Parse the arguments %%
s = parse_args(varargin);
if(isempty(s)) return; end

% Get the output extension
ext = getenv('FSF_OUTPUT_FORMAT');
if(~isempty(ext)) 
  s.UseMRIread = 1;
else
  ext = 'bhdr'; 
end
fprintf('UseMRIread = %d, ext = %s\n',s.UseMRIread,ext);

s = check_params(s);
if(isempty(s)) return; end

% This may be needed 
%s.PreStimWin = s.TER*floor(s.PreStimWin/s.TER);
%s.TotWin = s.TER*round(s.TotWin/s.TER);
outvolpath = fast_dirname(deblank(s.hvol));

sxa_print_struct(s,1);

TR  = s.TR;
TER = s.TER;
TW  = s.TotWin;
TPS = s.PreStimWin;
RmBaseline = s.RmBaseline;
RmTrend    = s.RmTrend;
QTrendFit  = s.QTrendFit;
RescaleTarget = s.RescaleTarget;
GammaFit = s.GammaFit;
gfDelta = s.gfDelta;
gfTau   = s.gfTau;
gfAlpha = s.gfAlpha;
nskip = s.nSkip;
firstslice = s.firstslice;
nslices = s.nslices;
TimeOffset = s.TimeOffset;
HanRadius = s.HanRad;
AcqOrder = s.AcqOrder;
SynthSeed = s.SynthSeed;


parfilelist = s.parlist;

instemlist = s.invollist;

tpxlist = s.tpxlist;

hstem     = s.hvol;
eresdir   = s.eresdir;
sigestdir = s.sigestdir;
pctstem   = s.pscvol;
fomnibusstem = s.fomnibusvol;
pomnibusstem = s.pomnibusvol;

if(~isempty(s.maskid))
  fprintf('INFO: loading mask %s\n',s.maskid);
  mask = MRIread(s.maskid);
  if(isempty(mask))
    fprintf('ERROR: could not load %s\n',s.maskid);
    return;
  end
  mask = mask.vol;
  nmasktot = length(find(mask));
  fprintf('INFO: mask has %d points\n',nmasktot);
else
  mask = [];
end


%-------------------------------------------------%
lastslice = firstslice + nslices - 1;
nruns = size(parfilelist,1);
Nfir = round(TW/TER);

if(SynthSeed < 0) SynthSeed = sum(100*clock); end
fprintf('SynthSeed = %10d\n',SynthSeed);

% Load Whitening matrix %
if(s.SecondPass)
   WhtnMtx = fmri_ldbfile(s.WhtnMtxFile);
   WhtnMtxSq = WhtnMtx'*WhtnMtx; %'
else
   WhtnMtx = [];
   WhtnMtxSq = [];
end

%-- Determine Number of Conditions across all runs----------%
% Get list of unqiue condition IDs for all runs %
condlist = [];
for run = 1:nruns
  par = fmri_ldpar(deblank(parfilelist(run,:)));
  if(isempty(par))
    fprintf('ERROR: reading parfile %s\n',deblank(parfilelist(run,:)));
    return;
  end
  condid = par(:,2);
  clrun = unique(par(:,2));
  condlist = unique([clrun; condlist]);
end

% Count number per condition

% Remove -1 and 0 %
ind = find(condlist ~= -1 & condlist ~= 0);
condlist = condlist(ind);

Nnnc = length(condlist); % excludes null %
Nc = Nnnc + 1;

fprintf(1,'Conditions Found (%d): ',Nnnc);
fprintf(1,'%2d ',condlist);
fprintf(1,'\n');

if(max(abs(diff(condlist))) > 1)
  fprintf('ERROR: the condition identifiers as found the the paradigm files\n');
  fprintf('       do not appear to be contiguous.\n');
  return;
end


% Check for holes in the list of condition numbers %
if(~isempty(find(diff(condlist)~=1)))
  fprintf(2,'fast_selxavg2: missing conditions\n');
  return;
end

% Count the number per condition %
parall = [];
Npercond= zeros(Nc,1);
for run = 1:nruns
  par = fmri_ldpar(deblank(parfilelist(run,:)));
  npc = [0];
  fprintf(1,'Run %2d: ',run);
  for c = condlist'  %'
    nc = length(find(par(:,2)==c)); 
    npc = [npc nc];
    fprintf(1,'%3d ',nc);
  end
  fprintf(1,'\n');
  Npercond = Npercond + npc'; %'
  parall = [parall; par];
end

% Setup Spatial Smoothing Kernal %
if(HanRadius > 1)
  HanFilter = fmri_hankernel(HanRadius);
end

SubSampRate = round(TR/TER);

% Get basic info from the first run %
instem = deblank(instemlist(1,:));
mri = MRIread(instem,1);
mristruct = mri.bhdr;
nrows = mri.volsize(1);
ncols = mri.volsize(2);
nslices = mri.volsize(3);
ntrs = mri.nframes;

%-----------------------------------------------------------------%
%--------------- Beginning of Slice Loop -------------------------%
%-----------------------------------------------------------------%
hmri = mri;
hoffsetmri = mri;
Fmri = mri;
pmri = mri;
rstd = mri;
rstd.vol = zeros(rstd.volsize);
ar1 = mri;
ar1.vol = zeros([ar1.volsize nruns]);

SumESSMtxRun = zeros(ntrs,ntrs,nruns);
NBrainVoxsRun = zeros(nruns);
SumESSMtx = 0;
NBrainVoxsTot = 0;
tic;
for slice = firstslice:lastslice
  fprintf(1,'Slice %d, %g --------------\n',slice,toc);
  SumXtX = 0;
  SumXtWX = 0;
  SumXtWY = 0;
  eres_ss = 0;
  Xfinal = [];

  for Pass = 1:2

    if(slice == firstslice | s.debug)
      if(Pass == 1)
        fprintf(1,'  First Pass (Accumulation), %g \n',toc);
      else
        fprintf(1,'  Second Pass (Residual Error Estimation), %g \n',toc);
      end
    end

    randn('state',SynthSeed); 

    if(Pass == 2)

      c = cond(SumXtWX);
      if(slice == firstslice | s.debug)
        fprintf('    Paradigm Condition: %g\n',c);
      end

      if(c > 10000000)
	fprintf('\n');
        fprintf('ERROR: paradigm is ill-conditioned (%g).\n',c);
        fprintf('Check your paradigm file for file for periodicities\n');
        fprintf('or for some event types that all ways follow other\n');
        fprintf('event types (or itself).\n');
        fname = sprintf('%s.sumxtx.mat',hstem);
	fprintf('Saving Xfinal to %s/%s\n',pwd,fname);
	save(fname,'SumXtWX','Xfinal');
	if(TR ~= TER)
	  fprintf('\n');
	  fprintf('It could also be due to the fact that the TER\n');
	  fprintf('does not equal the TR. SubTR estimation only works\n');
	  fprintf('when the event schedule has been optimized with\n');
	  fprintf('SubTR estimation in mind.\n');
	  fprintf('\n');
	end
	if(~isempty(tpxlist))
	  fprintf('\n');
	  fprintf('It could also be due to excluded time points.\n');
	  fprintf('Try running without excluded time points to .\n');
	  fprintf('see if the problem goes away.\n');
	  fprintf('\n');
	end
        return;
      end

      if(slice == firstslice | s.debug)
        indtask = 1:size(Xpar,2);
        invSumXtX = inv(SumXtWX);
        d = diag(invSumXtX);
        minvr  = min(1./d(indtask));
        meanvr = mean(1./d(indtask));
        maxvr  = max(1./d(indtask));
        fprintf('    Var Reduction Range: %g %g %g\n',minvr,meanvr,maxvr);
      end

      if(slice == firstslice | s.debug)
        fprintf(1,'     Computing hhat\n');
      end
      hCovMtx = inv(SumXtWX);
      hhat = hCovMtx*SumXtWY;
    end

    %-------------------- Run Loop ------------------------%
    Xfirall = [];
    ntrstot = 0;
    ntpxtot = 0;
    for run = 1:nruns

      if(slice == firstslice | s.debug)
        fprintf(1,'     Run %d/%d, %g \n',run,nruns,toc);
      end

      instem = deblank(instemlist(run,:));
      instemdir = fast_dirname(instem);

      % Get number of TRs in this run %
      if(~s.UseMRIread)
	[nslices nrows ncols ntrs] = fmri_bvoldim(instem);
      else
	% This might be slow for compressed nifti or mgz
	mri = MRIread(instem,1);
	if(~isempty(mri)) ntrs = mri.nframes;
	else              ntrs = 0;
	end
      end
      if(ntrs == 0)
	fprintf('ERROR: with %s\n',instem);
	return;
      end
      ntrstot = ntrstot + ntrs;
      nvoxs = nrows*ncols;

      % Time Point Exclusion %
      if(~isempty(tpxlist))
        TPExcludeFile = deblank(tpxlist(run,:));
        if(strcmp(TPExcludeFile,'noexcl')) TPExcludeFile = []; end
      else
        TPExcludeFile = [];
      end
      [indTPExcl indTPIncl] = fast_ldtpexcl(TPExcludeFile,TR,ntrs,nskip);
      ntpx = length(indTPExcl);
      if(slice == firstslice | s.debug)
        fprintf(1,'       Excluding %d Points: ',ntpx);
        fprintf(1,'%d ',indTPExcl);
        fprintf(1,'\n');
      end
      ntpxtot = ntpxtot + ntpx; 
      FF = min(indTPIncl); % First Frame

      Xdrift = [];
      if(s.PFOrder < 0)
        % Create Baseline/Trend Components of Convolution Matrix %
        Xbaseline = []; Xtrend = []; Xqtrend  = [];
        if(RmBaseline) Xbaseline = fast_baselinemtx(run,ntrs,nruns); end
        if(RmTrend)    Xtrend    = fast_trendmtx(run,ntrs,nruns); end
        if(QTrendFit)  Xqtrend   = fast_quadtrendmtx(run,ntrs,nruns); end
        Xdrift = [Xbaseline Xtrend Xqtrend];
      else
        Xdrift  = fast_polytrendmtx(run,ntrs,nruns,s.PFOrder);
      end

      if(s.nyqreg)
	Xnyq = ones(ntrs,1);
	Xnyq(2:2:end) = -1;
	Xnyq = Xnyq - mean(Xnyq); %Make sure it's demeaned
      else
	Xnyq = [];
      end
      
      % Load paradigm for this run %
      par = fmri_ldpar(deblank(parfilelist(run,:)));
      if(s.autostimdur)
	par = parboxcar(par,TER,[],ntrs*TR);
	if(isempty(par))
	  fprintf('ERROR: applying auto stimulus duration\n');
	  return;
	end
      end
      if(~isempty(s.stimdur))
	par = parboxcar(par,TER,s.stimdur);
	if(isempty(par))
	  fprintf('ERROR: applying stimulus duration\n');
	  return;
	end
      end
      
      % Compute Offset for Slice Timing %
      if(~isempty(AcqOrder))
        SliceDelay = fast_slicedelay(TR,nslices,slice,AcqOrder);
      else
        SliceDelay = 0;
      end

      % Adjust for Time Offset %
      par(:,1) = par(:,1) + TimeOffset;

      % Convert paradigm to FIR stimulus convolution matrix %
      Xfir = fmri_par2scm(par,Nc,SubSampRate*ntrs,TER,Nfir,TPS);

      % For Sub-TR Estimation %
      if(TR ~= TER)
        Xfirtmp = Xfir;
        nn = [1:SubSampRate:size(Xfirtmp,1)];
        Xfir = Xfirtmp(nn,:);
      end
      Xfirall = [Xfirall; Xfir];

      % Tranform for Fitting to Gamma Function or SPM HRF %
      if(GammaFit > 0)
        Xpar = fmri_scm2gcm(Xfir,Nnnc,TER,TPS,gfDelta,gfTau,gfAlpha);
        Navgs_per_cond = length(gfDelta);
      elseif(s.spmhrf > -1)
	tspmhrf = TER*[0:Nfir-1]'-TPS;
	hspmhrf = fast_spmhrf(tspmhrf);
	Aspmhrf = hspmhrf;
	dhspmhrf = hspmhrf;
	for nderiv = 1:s.spmhrf
	  % Divide by TER for gradient.
	  % Multiply by 2.6 to bring 1st deriv to amp of 1
	  dhspmhrf = 2.6*gradient(dhspmhrf)/TER;
	  Aspmhrf = [Aspmhrf dhspmhrf];
	end
	A = [];
	for c = 1:Nnnc
	  A = fast_blockdiag2(A,Aspmhrf);
	end
	Xpar = Xfir*A;
	Navgs_per_cond = s.spmhrf+1;
      else
        Xpar = Xfir;
        Navgs_per_cond = Nfir;
      end

      % Number of averages excluding the offsets and trends
      NTaskAvgs = Nnnc*Navgs_per_cond;

      if(~isempty(s.extreglist))
        extregstem = deblank(s.extreglist(run,:));
	if(exist(extregstem,'file'))
	  extreg = load(extregstem);
	else
	  extreg = fmri_ldbvolume(extregstem);
	  if(isempty(extreg))
	    fprintf('ERROR: could not load %s\n',extregstem);
	    return;
	  end
        end
        if(size(extreg,3)~=1) extreg = squeeze(extreg)'; %'
        else                  extreg = squeeze(extreg);
        end
        if(slice == firstslice & s.nextreg < 0)s.nextreg = size(extreg,2); end
        if(s.nextreg > size(extreg,2))
          fprintf('ERROR: %s does not have enough regressors\n',extregstem);
          return;
        end
        % Remove mean of External Regressor %
        extreg = extreg(:,1:s.nextreg);
        extreg = extreg - repmat(mean(extreg), [ntrs 1]);
        extreg = extreg./repmat(std(extreg), [ntrs 1]);
        if(s.extregorthog)
          extreg = ( eye(ntrs) - Xpar*inv(Xpar'*Xpar)*Xpar') * extreg;
        end
        z = zeros(size(extreg));
        extregrun = [repmat(z,[1 (run-1)]) extreg repmat(z,[1 (nruns-run)])];
      else
        extregrun = [];
      end

      % Create final Convolution Matrix for ith run %
      Xi = [Xpar Xdrift extregrun Xnyq];

      % Load or synthsize data %
      if(SynthSeed == 0)
        % Load the data for this slice %
	if(~s.UseMRIread)
	  [nrows ncols ntp fs ns endian bext] = fmri_bfiledim(instem);
	else
	  if(slice == firstslice)
	    mri = MRIread(instem,1);
	  end
	  nrows = mri.volsize(1);
	  ncols = mri.volsize(2);
	  ntp = mri.nframes;
	  ns = mri.volsize(3);
	end
	if(~isempty(s.TauMaxWhiten) & s.TauMaxWhiten ~= 0)
	  if(run == 1) ntprun = ntp;
	  else
	    if(ntp ~= ntprun)
	      fprintf(['ERROR: when using temporal whitening, all ' ...
		       'runs must have the same number of time ' ...
		       'points. Run 1 has %d time points whereas ' ...
		       'Run %d has %d time points.\n'],ntprun,run,ntp);
	      return;
	    end
	  end
	end
	if(~s.UseMRIread)
	  fname = sprintf('%s_%03d.%s',instem,slice,bext);
	  y = fmri_ldbfile(fname);
	else
	  if(slice == 0 & Pass == 1)  
	    fprintf('       Loading %s\n',instem);
	    ymri(run) = MRIread(instem); 
	  end
	  y = ymri(run).vol(:,:,slice+1,:);
	  y = permute(y,[1 2 4 3]); % instead of squeeze
	  %y = squeeze(ymri(run).vol(:,:,slice+1,:));
	end
      else
        %fprintf(1,'       Synthesizing Data for Slice %d \n',slice);
        y = randn(nrows, ncols, ntrs);
      end

      % Exlude Points %
      Xi(indTPExcl,:) = 0;
      y(:,:,indTPExcl)  = 0;

      % Spatial Smoothing with In-Plane Hanning Window %
      if(HanRadius > 1)
        if(slice == firstslice | s.debug)
          fprintf(1,'       Spatial Smoothing, HanRad = %g\n',HanRadius);
        end
        y = fmri_spatfilter(y,HanFilter);
      end

      % Reshape to a more convenient form %
      y = reshape(y, [nrows*ncols ntrs])'; %'

      % Global rescale of functional data %
      if(RescaleTarget > 0)
        %MeanValFile = sprintf('%s.meanval',instem);
        MeanValFile = sprintf('%s/global.meanval.dat',instemdir);
        [RescaleFactor MeanVal]=fast_rescalefactor(MeanValFile, RescaleTarget);
        %fprintf(1,'       Rescaling Global Mean %g,%g,%g\n',...
     	%        MeanVal,RescaleTarget,RescaleFactor);
        y = RescaleFactor * y;
	%fprintf('MeanVal = %g\n',MeanVal);
      else
        RescaleFactor = 1;
      end
      %fprintf('RescaleFactor = %g\n',RescaleFactor);
      %fprintf('RescaleTarget = %g\n',s.RescaleTarget);

      if(s.loginput) 
        fprintf('INFO: computing log of input\n');
        y = log(abs(y)+2); 
      end

      % Load per-run whitening matrix here, mult by Xi and y
      if(s.WhitenFlag)
        fname = deblank(s.WhtnMtxFile(run,:));
        WW = load(fname);
        WhtnMtx = WW.W;
        WhtnMtxSq = WhtnMtx'*WhtnMtx; %'
      else
        % Whitening Filter %
        if(isempty(s.WhtnMtxFile)) 
          WhtnMtx   = eye(ntrs); 
          WhtnMtxSq = eye(ntrs); 
        else
          [wr wc] = size(WhtnMtx);
          if(wr ~= ntrs)
            fprintf(2,'ERROR: Whitening Matrix is %d x %d, ntrs = %d\n',...
                    wr,wc,ntrs);
            return;
          end
        end
      end

      if(Pass == 1)
        % Accumulate XtX and XtWY %
        SumXtX  = SumXtX  + Xi'*Xi; %'
        SumXtWX = SumXtWX + Xi'*WhtnMtxSq*Xi; %'
        SumXtWY = SumXtWY + (Xi'*WhtnMtxSq)*y;  %'
        Xfinal = [Xfinal; Xi];
      end

      if(Pass == 2)
        sigest = Xi*hhat;  % Compute Signal Estimate %

        % Compute Residual Error %
        if(isempty(s.WhtnMtxFile)) 
          eres = y - sigest; % Unwhitened
        else
          eres = WhtnMtx*(y - sigest); % Whitened
        end

        % Accumulate Sum of Squares of  Residual Error %
        eres_ss = eres_ss + sum(eres.^2);
        
        if(~isempty(sigestdir))
          % Save Signal Estimate (Partial Model Fit) %
  	  pmf = Xi(:,1:NTaskAvgs)*hhat(1:NTaskAvgs,:);
          fname = sprintf('%s/s%03d_%03d.bfloat',sigestdir,run,slice);
          tmp = reshape(pmf', [nrows ncols ntrs])/RescaleFactor; %'
  	  fmri_svbfile(tmp,fname);
          bhdrfile = sprintf('%s/s%03d.bhdr',sigestdir,run);
	  fast_svbhdr(mristruct, bhdrfile, 0);
  	  clear tmp pmf;
        end

        if(~isempty(eresdir))
          % Save (Whitened) Residual Error %
          fname = sprintf('%s/e%03d_%03d.bfloat',eresdir,run,slice);
          tmp = reshape(eres', [nrows ncols ntrs])/RescaleFactor; %'
  	  fmri_svbfile(tmp,fname);
          bhdrfile = sprintf('%s/e%03d.bhdr',eresdir,run);
	  fast_svbhdr(mristruct, bhdrfile, 0);
        end

        if(~isempty(s.acfdir))
          % Compute and save ACF
          fname = sprintf('%s/acf%03d_%03d.bfloat',s.acfdir,run,slice);
	  if(slice == firstslice | s.debug)
	    fprintf('INFO: computing acf %g\n',toc);
	  end
	  acf = fast_acorr(eres);
          tmp = reshape(acf', [nrows ncols ntrs]);
  	  fmri_svbfile(tmp,fname);
        end

	if(s.NewWhitenFlag)
	  fprintf('  Run %d, Computing ar1\n',run);
	  indBrain = find(mask(:,:,slice+1));
	  if(~isempty(indBrain))
	    ar1mask = sum(eres(1:end-1,indBrain).*eres(2:end,indBrain))./...
		sum(eres(:,indBrain).^2);
	    ar1slice = zeros([ncols nrows]);
	    ar1slice(indBrain) = ar1mask;
	    ar1.vol(:,:,slice+1,run) = ar1slice;
	  end
	end
	
        if(~isempty(s.ErrCovMtxStem))
 	  if(s.SegBrainAir)
	    if(isempty(mask))
              MeanValFile = sprintf('%s.meanval',instem);
              [tmp MeanVal] = fast_rescalefactor(MeanValFile,101);
   	      indBrain = find( y(FF,:) > .75*MeanVal );
            else
              indBrain = find(mask(:,:,slice+1));
            end
          else
  	    indBrain = [1 nvoxs];
          end
          NBrainVoxs = length(indBrain);
          if(run==1) fprintf(1,' NBrainVoxs = %d\n',NBrainVoxs);end

	  if(NBrainVoxs > 0)
    	    %fprintf(1,'       Computing Err SS Matrix\n');
            ESSMtx = eres(:,indBrain) * eres(:,indBrain)'; % '
            SumESSMtx = SumESSMtx + ESSMtx;
            NBrainVoxsTot = NBrainVoxsTot + NBrainVoxs;
            SumESSMtxRun(:,:,run) = SumESSMtxRun(:,:,run) + ESSMtx;
            NBrainVoxsRun(run) = NBrainVoxsRun(run) + NBrainVoxs;
          end
        end % isempty(ECVMStem)

      end % if(Pass == 2)

    end % Loop over runs %
  end % Loop over Pass %

  
  
  % Total Number of Averages Computed %
  Navgs_tot = size(SumXtX,1);

  % Residual Error Forming Matrix 
  R = eye(size(Xfinal,1)) - Xfinal*inv(Xfinal'*Xfinal)*Xfinal'; 

  % Total Degrees of Freedom
  DOF = trace(R);

  %fprintf(1,'  Computing Residual Error Std\n');
  eres_var = eres_ss/DOF;
  eres_std = sqrt(eres_var);

  % -------- Convert to selavg format -------------- %
  hhattmp = hhat(1:NTaskAvgs,:); %Remove offset and baseline 
  hhattmp = [zeros(Navgs_per_cond,nvoxs); hhattmp]; % Add zero for cond 0
  hhattmp2 = reshape(hhattmp,[Navgs_per_cond Nc nvoxs]);

  hstd = sqrt( (diag(hCovMtx).*diag(SumXtX)) * eres_std.^2);
  hstdtmp = hstd(1:NTaskAvgs,:); % Remove offset and baseline
  hstdtmp = [repmat(eres_std, [Navgs_per_cond 1]); hstdtmp]; % Add 0 for cond 0
  hstdtmp2 = reshape(hstdtmp,[Navgs_per_cond Nc nvoxs]);

  %--- Merge Averages and StdDevs ---%
  tmp = zeros(Navgs_per_cond,2,Nc,nvoxs);
  tmp(:,1,:,:) = hhattmp2;
  tmp(:,2,:,:) = hstdtmp2;
  tmp = reshape(tmp,[Navgs_per_cond*2*Nc nrows ncols ]);
  tmp = permute(tmp,[2 3 1]);

  % Save in selxavg format %
  if(~s.UseMRIread)
    fname = sprintf('%s_%03d.bfloat',hstem,slice);
    if(slice == firstslice | s.debug)
      fprintf(1,'  Saving data to %s \n',fname);
    end
    fmri_svbfile(tmp,fname);
    %if(s.UseMRIread) fast_svbhdr(mristruct,hstem,1); end
  else
    hmri.vol(:,:,slice+1,:) = tmp;
  end
    
  % Save the mean image %
  if(RmBaseline)
    hoffset = reshape(hhat(NTaskAvgs+1,:), [nrows ncols ]); % From 1st Run
    indz = find(hoffset==0);
    hoffset(indz) = 1;
    fname = sprintf('%s-offset_%03d.bfloat',hstem,slice);
    if(slice == firstslice | s.debug)
      fprintf(1,'  Saving offset to %s \n',fname);
    end
    if(~s.UseMRIread)
      fmri_svbfile(hoffset,fname);
      stem = sprintf('%s-offset',hstem);
      %if(s.UseMRIread) fast_svbhdr(mristruct,stem,1); end
    else
      hoffsetmri.vol(:,:,slice+1) = hoffset;
    end
  end

  % Save Percent Signal Chanage %
  if(~isempty(pctstem))
    fprintf('Saving pct\n');
    hofftmp = hoffset;
    ind = find(hofftmp==0);
    hofftmp(ind) = 1e10;
    tmp = 100 * tmp ./ repmat(hofftmp, [1 1 size(tmp,3)]);
    fname = sprintf('%s_%03d.bfloat',pctstem,slice);
    fprintf(1,'  Saving Percent Signal Change to %s \n',fname);
    fmri_svbfile(tmp,fname);
    fname = sprintf('%s-offset_%03d.bfloat',pctstem,slice);
    fprintf(1,'  Saving pctsigch offset to %s \n',fname);
    fmri_svbfile(ones(nrows,ncols),fname);
  end

  % Save in beta format %
  if(~isempty(s.betavol))
    tmp = hhat;
    ntmp = size(tmp,1);
    tmp = reshape(tmp',[nrows ncols ntmp]); %';
    err = fast_svbslice(tmp,s.betavol,slice,'',mristruct);
    if(err) 
      fprintf('ERROR: saving %s\n',s.betavol);
      return;
    end
    clear tmp;
  end
    
  tmp = eres_var;
  tmp = reshape(tmp',[nrows ncols]);
  rstd.vol(:,:,slice+1) = sqrt(tmp);

  % Omnibus Significance Test %
  if(~isempty(fomnibusstem) | ~isempty(pomnibusstem))

    R = eye(NTaskAvgs,Navgs_tot);
    q = R*hhat;
    if(size(q,1) > 1)  qsum = sum(q); % To get a sign %
    else               qsum = q;
    end

    if(NTaskAvgs == 1)
      Fnum = inv(R*hCovMtx*R') * (q.^2) ;  %'
    else
      Fnum = sum(q .* (inv(R*hCovMtx*R') * q));  %'
    end
    Fden = NTaskAvgs*(eres_std.^2);
    ind = find(Fden == 0);
    Fden(ind) = 10^10;
    F = sign(qsum) .* Fnum ./ Fden;
    if(~isempty(fomnibusstem))
      fname = sprintf('%s_%03d.bfloat',fomnibusstem,slice);
      tmp = reshape(F,[nrows ncols]);
      if(~s.UseMRIread) 
	fmri_svbfile(tmp,fname);
	%if(s.UseMRIread) fast_svbhdr(mristruct,fomnibusstem,1); end
      else
	Fmri.vol(:,:,slice+1) = tmp;
      end
    end
    if(~isempty(pomnibusstem))
      if(slice == firstslice | s.debug)
        fprintf('INFO: performing FTest on omnibus\n');
        fprintf('      NOTE: if this hangs, try running selxavg-sess\n');
        fprintf('      with the -noomnibus flag.\n');
	if(s.UseMRIread)
	  fprintf('Using MRIread(), which can slow things down.\n');
	end
      end
      p = sign(F) .* FTest(NTaskAvgs,DOF,abs(F));
      indz = find(p==0);
      p(indz) = 1;
      p = sign(p).*(-log10(abs(p)));
      tmp = reshape(p,[nrows ncols]);
      if(~s.UseMRIread) 
	fname = sprintf('%s_%03d.bfloat',pomnibusstem,slice);
	fmri_svbfile(tmp,fname);
	%if(s.UseMRIread) fast_svbhdr(mristruct,pomnibusstem,1); end
      else
	pmri.vol(:,:,slice+1) = tmp;
      end
    end
  end

  if(~isempty(s.snrdir))
    p = zeros(Navgs_tot,1);
    p(1:NTaskAvgs) = 1;
    P = diag(p);
    sigvar = sum(hhat .* ((P'*SumXtX*P) * hhat))/(ntrstot);  %'

    sigvar = reshape(sigvar,[nrows ncols]);
    fname = sprintf('%s/sigvar_%03d.bfloat',s.snrdir,slice);
    fmri_svbfile(sigvar,fname);
    fname = sprintf('%s/sigstd_%03d.bfloat',s.snrdir,slice);
    fmri_svbfile(sqrt(sigvar),fname);
    
    resvar = reshape(eres_std.^2,[nrows ncols ]);
    fname = sprintf('%s/resvar_%03d.bfloat',s.snrdir,slice);
    fmri_svbfile(resvar,fname);
    fname = sprintf('%s/resstd_%03d.bfloat',s.snrdir,slice);
    fmri_svbfile(sqrt(resvar),fname);
    fname = sprintf('%s/resstdpct_%03d.bfloat',s.snrdir,slice);
    fmri_svbfile(100*sqrt(resvar)./hoffset,fname);

    ind0 = find(resvar==0);
    resvar(ind0) = 1;
    snrvar = sigvar./resvar;
    snrvar(ind0) = 0;
    fname = sprintf('%s/snrvar_%03d.bfloat',s.snrdir,slice);
    fmri_svbfile(snrvar,fname);
    fname = sprintf('%s/snrstd_%03d.bfloat',s.snrdir,slice);
    fmri_svbfile(sqrt(snrvar),fname);
    
  end

end % Loop over slices 
%------------------------------------------------------------%

if(s.UseMRIread)
  fname = sprintf('%s.%s',hstem,ext);
  MRIwrite(hmri,fname);
  fname = sprintf('%s-offset.%s',hstem,ext);
  MRIwrite(hoffsetmri,fname);
  fname = sprintf('%s.%s',pomnibusstem,ext);
  MRIwrite(pmri,fname);
  fname = sprintf('%s.%s',fomnibusstem,ext);
  MRIwrite(Fmri,fname);
end

if(s.NewWhitenFlag)
  % Load the offest back in
  fname = sprintf('%s-offset',hstem);
  hoffset = MRIread(fname);

  outvolpath = fast_dirname(deblank(s.hvol));
  fname = sprintf('%s/ar1.mgh',outvolpath);
  MRIwrite(ar1,fname);
  ar1mn = mean(ar1.vol,4);

  acfseg = mri;
  acfseg.vol = fast_acfseg(hoffset.vol,10,mask);
  fname = sprintf('%s/acfseg.nii',outvolpath);
  MRIwrite(acfseg,fname);
  segunique = unique(sort(acfseg.vol(:)));
  nunique = length(segunique);
  for nth = 0:nunique-1
    indk = find(acfseg.vol == nth);
    ar1seg(nth+1) = mean(ar1mn(indk));
  end
  acffile = sprintf('%s/acf.mat',outvolpath);
  save(acffile,'ar1seg','acfseg','ar1');
end

xfile = sprintf('%s/X.mat',outvolpath);
pfOrder = s.PFOrder;
nExtReg = 0; if(s.nextreg > 0) nExtReg = s.nextreg; end
tPreStim = s.PreStimWin;
TimeWindow = TW;
nyqreg = s.nyqreg;
spmhrf = s.spmhrf;
autostimdur = s.autostimdur;
stimdur = s.stimdur;
fprintf('INFO: saving meta to %s\n',xfile);
save(xfile,'Xfinal','Nnnc','pfOrder','nExtReg',...
     'nruns','Navgs_per_cond','TimeWindow','tPreStim','TR','TER',...
     'gfDelta','gfTau','gfAlpha','tpxlist','RescaleFactor',...
     'RescaleTarget','nyqreg','spmhrf','autostimdur','stimdur','-v4');

%-- Save ECovMtx for each run individually --%
if(s.SaveErrCovMtx) 
  for run = 1:nruns,
    fname = sprintf('%s-ecvm_%03d.bfloat',s.hvol,run-1);
    tmp = SumESSMtxRun(:,:,run)/NBrainVoxsRun(run);
    fmri_svbfile(tmp,fname);
  end
end

% Save ErrCovMtxs across all runs %/
if(~isempty(s.ErrCovMtxStem))
  fprintf(1,'NBrainVoxsTot: %d\n',NBrainVoxsTot);
  fname = sprintf('%s.bfloat',s.ErrCovMtxStem);
  ErrCovMtx = SumESSMtx/NBrainVoxsTot;
  fmri_svbfile(ErrCovMtx,fname);
end

fname = sprintf('%s/rstd.%s',outvolpath,ext);
MRIwrite(rstd,fname);
rvar = rstd;
rvar.vol = rstd.vol.^2;
fname = sprintf('%s/rvar.%s',outvolpath,ext);
MRIwrite(rvar,fname);

% Save the .dat file %
fname = sprintf('%s.dat',hstem);
SumXtXTmp  = SumXtX( 1:NTaskAvgs, 1:NTaskAvgs);
hCovMtxTmp = hCovMtx(1:NTaskAvgs, 1:NTaskAvgs);
hd = fmri_hdrdatstruct;
hd.TR  = TR;
hd.TER = TER;
hd.TimeWindow = TW;
hd.TPreStim = TPS;
hd.Nc = Nc;
hd.Nh = Navgs_per_cond;
hd.Nnnc = Nnnc;
hd.DOF= DOF;
hd.Npercond= Npercond;
hd.Nruns = nruns;
hd.Ntp = ntrstot;
hd.Nrows = nrows;
hd.Ncols = ncols;
hd.Nskip = nskip;
if(s.PFOrder < 0)
  hd.DTOrder = RmBaseline+RmTrend+QTrendFit;
else
  hd.DTOrder = s.PFOrder + 1;
end
hd.RescaleFactor = RescaleFactor;
hd.HanningRadius = 0.0;
hd.BrainAirSeg = 0;
hd.GammaFit = GammaFit ;
hd.gfDelta  = gfDelta;
hd.gfTau    = gfTau;
if(s.spmhrf > -1) % Hack
  hd.GammaFit = s.spmhrf + 1;
  hd.gfDelta = -1*ones(1,s.spmhrf + 1);
  hd.gfTau   = -1*ones(1,s.spmhrf + 1);
end;
%hd.gfAlpha      = gfAlpha; % This is not saved yet
hd.NullCondId    = 0;
hd.SumXtX        = SumXtXTmp;
hd.nNoiseAC      = 0;
hd.CondIdMap     = [0:Nc-1];
hd.hCovMtx       = hCovMtxTmp;
hd.WhitenFlag        = s.WhitenFlag;
hd.runlist  = getrunlist(s.invollist);
hd.funcstem      = basename(deblank(s.invollist(1,:)));
hd.parname       = s.parname;
if(~isempty(s.extreglist))
  hd.extregstem  = basename(deblank(s.extreglist(1,:)));
  hd.nextreg  = s.nextreg;
  hd.extortho = s.extregorthog;
end

fmri_svdat3(fname,hd);

if(~isempty(pctstem))
  fname = sprintf('%s.dat',pctstem);
  fmri_svdat3(fname,hd);
end

%------------------------------------------------------------%
if(s.AutoWhiten)

  fprintf('Computing Whitening Matrix\n');
  % Note that this can be done with the following function:
  %W2 = fast_ecvm2wmtx(ErrCovMtx);
  % Also note that TauMax and PctWhiten are not used:).
  % Also note that this produces a non-causal whitening matrix.
  
  s.NMaxWhiten = round(s.TauMaxWhiten/s.TR);
  fprintf('TauMax = %g, NMax = %d, Pct = %g\n',s.TauMaxWhiten,...
          s.NMaxWhiten,s.PctWhiten);
  nf = size(ErrCovMtx,1);

  if(1)
    AutoCor = fast_cvm2acor(ErrCovMtx,1);
    [cnd mineig S] = fast_acfcond(AutoCor);
    ncnd = 1;
    while(mineig < 0 | cnd > 100)
      fprintf('ncnd = %d, mineig = %g, cnd = %g\n',ncnd,mineig,cnd);
      AutoCor = AutoCor .* tukeytaper(nf);
      [cnd mineig S] = fast_acfcond(AutoCor);
      ncnd = ncnd + 1;
    end
    fprintf('ncnd = %d, mineig = %g, cnd = %g\n',ncnd,mineig,cnd);
    [utmp stmp vtmp] = svd(S);
    W = utmp * inv(sqrt(stmp)) * vtmp';
  end

  whtnmtxfile = sprintf('%s-whtnmtx.bfloat',s.hvol);
  fmri_svbfile(W,whtnmtxfile);
  autocorfile = sprintf('%s-autocor.dat',s.hvol);
  fid = fopen(autocorfile,'w');
  fprintf(fid,'%8.4f\n',AutoCor);
  fclose(fid);

  
  fprintf('Starting Whitening Stage\n');
  nvarargin = length(varargin);
  varargin{nvarargin+1} = '-noautowhiten';
  varargin{nvarargin+2} = '-whtnmtx';
  varargin{nvarargin+3} = whtnmtxfile;
  fprintf('\n\n\n\nStarting recursive call to fast_selxavg\n');
  r = fast_selxavg(varargin{:});
  fprintf('Recursive call to fast_selxavg finished\n\n\n');


  % This is the old way to compute the whitening matrix %
  % Condition Error Covariance Matrix %
  %[ErrCovMtxRgl neig pve] = fast_cvm_condrgl(ErrCovMtx,s.MinCond);
  % Normalze Regularized Error Covariance Matrix %
  %ErrCovMtxRglNorm = fast_cvm_normalize(ErrCovMtxRgl);
  %ErrCovMtxRglNorm = ErrCovMtxRgl/mean(diag(ErrCovMtxRgl));
  %[ErrCovMtxRglNorm niters alpha] = fast_cvm_normrgl(ErrCovMtxRgl,s.MinCond);
  %fprintf('Final Cond = %g, niters = %d, alpha = %g\n',...
  %        cond(ErrCovMtxRglNorm),niters,alpha);
  %fprintf('MinCond = %g, FinalCond = %g, NEig = %d, PVE = %g\n',...
  %         s.MinCond,cond(ErrCovMtxRglNorm),neig,pve);
  % Compute and save whitening matrix %
  %W = inv(ErrCovMtxRglNorm);


end
%------------------------------------------------------------%



fprintf(1,'Done %g\n',toc);

r = 0;

return;
%---\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%-----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%


%------------- Print Usage ---------------------%
function print_usage(dummy)

  fprintf(1,'USAGE:\n');
  fprintf(1,'  fast_selxavg2\n');
  fprintf(1,'     -i   invol ... \n');
  fprintf(1,'     -p   parfile ... \n');
  fprintf(1,'     -parname  parfile \n');
  fprintf(1,'     -extreg   extregfile \n');
  fprintf(1,'     -nextreg  number of external regressors to use\n');
  fprintf(1,'     -parname  parfile \n');
  fprintf(1,'     -tpx tpxfile ... \n');
  fprintf(1,'     -whtmtx whitening matrix file \n');
  fprintf(1,'     -o   hdrstem \n');
  fprintf(1,'     -psc pscstem \n');
  fprintf(1,'     -fomnibus stem \n');
  fprintf(1,'     -pomnibus stem \n');
  fprintf(1,'     -TR   TR\n');
  fprintf(1,'     -TER  TER\n');
  fprintf(1,'     -timewindow totwin  \n');
  fprintf(1,'     -prewindow  prewin  \n');
  fprintf(1,'     -nobaseline  \n');
  fprintf(1,'     -detrend  \n');
  fprintf(1,'     -qtrendfit  \n');
  fprintf(1,'     -rescale  target \n');
  fprintf(1,'     -nskip  n \n');
  fprintf(1,'     -hanrad radius \n');
  fprintf(1,'     -fwhm   width \n');
  fprintf(1,'     -ipr    inplaneres \n');
  fprintf(1,'     -gammafit delta tau \n');
  fprintf(1,'     -timeoffset t \n');
  fprintf(1,'     -acqorder  <linear or interleaved> \n');
  fprintf(1,'     -firstslice sliceno : 0 \n');
  fprintf(1,'     -nslices    nslices : auto \n');
  fprintf(1,'     -eresdir    dir \n');
  fprintf(1,'     -sigestdir  dir \n');
  fprintf(1,'     -synth      seed \n');
  fprintf(1,'     -cfg        file \n');

return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = sxa_struct
  s.invollist      = '';
  s.parlist        = '';
  s.autostimdur    = 0; % extract stim durations from par
  s.stimdur        = []; % stimulus duration
  s.nruns          = 0;
  s.parname        = '';
  s.extregfile     = '';
  s.extreglist     = '';
  s.nextreg        = -1;
  s.extregorthog   =  0;
  s.tpxlist        = '';
  s.WhtnMtxFile     = '';
  s.AutoWhiten     = 0;
  s.NoAutoWhiten   = 0;
  s.SecondPass     = 0;
  s.TauMaxWhiten   = 0;
  s.NMaxWhiten     = 0;
  s.PctWhiten      = 0;
  s.LPFFlag        = 0;
  s.HPF            = [];
  s.WhitenFlag     = 0;
  s.NewWhitenFlag  = 0;
  s.maskid         = []; % for whitening only
  s.hvol           = '';
  s.betavol        = '';
  s.fomnibusvol    = '';
  s.pomnibusvol    = '';
  s.ErrCovMtxStem  = '';
  s.SaveErrCovMtx  = 0;
  s.pscvol   = '';
  s.TR    = '';
  s.TER    = '';
  s.TotWin      = '';
  s.PreStimWin  = 0;
  s.PostStimWin = '';
  s.SegBrainAir = 1;
  s.RmBaseline = 1;
  s.RmTrend    = 0;
  s.QTrendFit  = 0;
  s.PFOrder    = -1;
  s.RescaleTarget = 0;
  s.nSkip  = 0;
  s.FWHM = 0;
  s.InPlaneRes = 0;
  s.HanRad = 0;
  s.gfDelta = [];
  s.gfTau = [];
  s.gfAlpha = 2;
  s.spmhrf = -1;
  s.TimeOffset = 0;
  s.AcqOrder = '';
  s.SynthSeed = 0;
  s.cfgfile = '';
  s.verbose = 0;
  s.firstslice = 0;
  s.nslices    = -1;
  s.eresdir    = '';
  s.sigestdir  = '';
  s.acfdir  = '';
  s.snrdir  = '';
  s.debug = 0;
  s.loginput = 0;
  s.funcstem = '';
  s.nyqreg = 0; % nyquist regressor
  s.UseMRIread = 0;
return;

%--------------------------------------------------%
% Parse the arguments from the config file %
function argscfg = parse_cfg(args)
  argscfg = args;
  cfgfile = '';
  nargs = length(args);
  narg = 1;
  while(narg <= nargs)
    flag = deblank(args{narg});
    narg = narg + 1;
    if(strcmp(flag,'-cfg'))
      arg1check(flag,narg,nargs);
      cfgfile = args{narg};
      break;
    end
  end

  if(~isempty(cfgfile))
    fid = fopen(cfgfile,'r');
    if(fid == -1)
      fprintf(2,'ERROR: cannot open %s\n',cfgfile);
      argscfg = []; return;
    end
    [s n] = fscanf(fid,'%s',1);
    while(n ~= 0)
      nargs = nargs + 1;;
      argscfg{nargs} = s;
      [s n] = fscanf(fid,'%s',1);
    end
    fclose(fid);
  end

return

%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = sxa_struct;
  inputargs = parse_cfg(varargin{1});
  ninputargs = length(inputargs);

  narg = 1;
  while(narg <= ninputargs)

    flag = deblank(inputargs{narg});
    narg = narg + 1;
    %fprintf(1,'Argument: %s\n',flag);
    if(~isstr(flag))
      flag
      fprintf(1,'ERROR: All Arguments must be a string\n');
      error;
    end

    switch(flag)

      case '-i',
        arg1check(flag,narg,ninputargs);
        s.invollist = strvcat(s.invollist,inputargs{narg});
        narg = narg + 1;

      case '-p',
        arg1check(flag,narg,ninputargs);
        s.parlist = strvcat(s.parlist,inputargs{narg});
        narg = narg + 1;

      case '-extreg',
       % Version5: this now has two args
       %arg1check(flag,narg,ninputargs);
        arg2check(flag,narg,ninputargs);
        s.extregfile = inputargs{narg};
        narg = narg + 1;
        s.nextreg = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;
     
     case '-extregorthog',
        s.extregorthog = 1;

      case '-nextreg',
        arg1check(flag,narg,ninputargs);
        s.nextreg = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-parname'},
        arg1check(flag,narg,ninputargs);
        s.parname = inputargs{narg};
        narg = narg + 1;

      case '-stimdur',
        arg1check(flag,narg,ninputargs);
        s.stimdur = sscanf(inputargs{narg},'%g',1);
        narg = narg + 1;

      case {'-autostimdur'},
        s.autostimdur = 1;

      case {'-noautostimdur'},
        s.autostimdur = 0;

      case {'-tpx','-tpexclfile'}
        arg1check(flag,narg,ninputargs);
        s.tpxlist = strvcat(s.tpxlist,inputargs{narg});
        narg = narg + 1;

      case {'-whtnmtx'}
        arg1check(flag,narg,ninputargs);
        s.WhtnMtxFile = strvcat(s.WhtnMtxFile,inputargs{narg});
        narg = narg + 1;

      case {'-autowhiten'} % Arg is minimum condition for regularization
        arg1check(flag,narg,ninputargs);
        s.TauMaxWhiten = sscanf(inputargs{narg},'%f',1);
        if(s.TauMaxWhiten > 0)
          s.AutoWhiten = 1;
          s.SaveErrCovMtx = 1;
        end
        narg = narg + 1;
     
     case {'-noautowhiten'} % To ease recursive calls
        s.NoAutoWhiten = 1;
        s.AutoWhiten = 0;
        s.TauMaxWhiten = 0;
        s.SecondPass = 1;

     case {'-acfseg'} % New Whiten
        s.NewWhitenFlag  = 1;
     
     case {'-mask'},
        arg1check(flag,narg,ninputargs);
        s.maskid = inputargs{narg};
        narg = narg + 1;

      case {'-o','-h'},
        arg1check(flag,narg,ninputargs);
        s.hvol = inputargs{narg};
        narg = narg + 1;

      case {'-beta'},
        arg1check(flag,narg,ninputargs);
        s.betavol = inputargs{narg};
        narg = narg + 1;

      case {'-fomnibus'},
        arg1check(flag,narg,ninputargs);
        s.fomnibusvol = inputargs{narg};
        narg = narg + 1;

      case {'-pomnibus'},
        arg1check(flag,narg,ninputargs);
        s.pomnibusvol = inputargs{narg};
        narg = narg + 1;

      case {'-ecovmtx'},
        arg1check(flag,narg,ninputargs);
        s.ErrCovMtxStem = inputargs{narg};
        narg = narg + 1;

      case {'-svecovmtx','-sverrcovmtx','-svecvm','-svacf'},
        s.SaveErrCovMtx = 1;

      case {'-psc','-percent'},
        arg1check(flag,narg,ninputargs);
        s.pscvol = inputargs{narg};
        narg = narg + 1;

      case {'-TR'}
        arg1check(flag,narg,ninputargs);
        s.TR = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-TER'}
        arg1check(flag,narg,ninputargs);
        s.TER = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-timewindow','-totwin','-tw'}
        arg1check(flag,narg,ninputargs);
        s.TotWin = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-prewindow','-prewin','-prestim'}
        arg1check(flag,narg,ninputargs);
        s.PreStimWin = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-postwindow','-postwin','-poststim'}
        arg1check(flag,narg,ninputargs);
        s.PostStimWin = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-timeoffset'}
        arg1check(flag,narg,ninputargs);
        s.TimeOffset = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-basegment'}
        s.SegBrainAir = 1;
  
      case {'-nobasegment'}
        s.SegBrainAir = 0;
  
      case {'-nobaseline'}
        s.RmBaseline = 0;
  
      case {'-baseline'}
        s.RmBaseline = 1;
  
      case {'-detrend'}
        s.RmTrend = 1;
  
      case {'-qtrendfit'}
        s.QTrendFit = 1;
  
      case {'-lpf'}
        s.LPFFlag = 1;
  
      case {'-hpf'}
        arg2check(flag,narg,ninputargs);
        s.HPF = sscanf(inputargs{narg},'%f %f',1);
        narg = narg + 1;

      case {'-polyfit'}
        arg1check(flag,narg,ninputargs);
        s.PFOrder = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-rescale'}
        arg1check(flag,narg,ninputargs);
        s.RescaleTarget = sscanf(inputargs{narg},'%f',1);
	fprintf('RescaleTarget = %g\n',s.RescaleTarget);
        narg = narg + 1;

      case {'-nskip'}
        arg1check(flag,narg,ninputargs);
        s.nSkip = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-hanrad'}
        arg1check(flag,narg,ninputargs);
        s.HanRad = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-fwhm'}
        arg1check(flag,narg,ninputargs);
        s.FWHM = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-ipr'}
        arg1check(flag,narg,ninputargs);
        s.InPlaneRes = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-gammafit'}
        arg2check(flag,narg,ninputargs);
        gfDelta = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;
        gfTau   = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;
        s.gfDelta = [s.gfDelta gfDelta];
        s.gfTau   = [s.gfTau   gfTau];

      case {'-gammaexp'}
        arg1check(flag,narg,ninputargs);
        s.gfAlpha = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case '-spmhrf', % Argument is number of derivatives
        arg1check(flag,narg,ninputargs);
        s.spmhrf = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-acqorder',
        arg1check(flag,narg,ninputargs);
        s.AcqOrder = inputargs{narg};
        narg = narg + 1;

      case {'-firstslice', '-fs'}
        arg1check(flag,narg,ninputargs);
        s.firstslice = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-nslices', '-ns'}
        arg1check(flag,narg,ninputargs);
        s.nslices = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-eresdir',
        arg1check(flag,narg,ninputargs);
        s.eresdir = inputargs{narg};
        narg = narg + 1;

      case '-acfdir',
        arg1check(flag,narg,ninputargs);
        s.acfdir = inputargs{narg};
        narg = narg + 1;

      case '-snrdir',
        arg1check(flag,narg,ninputargs);
        s.snrdir = inputargs{narg};
        narg = narg + 1;

      case {'-sigestdir','-signaldir'}
        arg1check(flag,narg,ninputargs);
        s.sigestdir = inputargs{narg};
        narg = narg + 1;

      case '-cfg',
        % This is actually handled by parse_cfg
        arg1check(flag,narg,ninputargs);
        narg = narg + 1;

      case '-synth', 
        arg1check(flag,narg,ninputargs);
        s.SynthSeed = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-nyqreg',
        s.nyqreg = 1;

      case '-verbose',
        s.verbose = 1;

      case {'-mriread'}
        s.UseMRIread = 1;
	fprintf('Using MRIread()\n');
  
      case '-log',
        s.loginput = 1;

      % ignore these guys %
      case {'-monly', '-nullcondid','-umask','-sveres','-svsignal','-svsnr','-acfbins'},
        arg1check(flag,narg,ninputargs);
        narg = narg + 1;

      case {'-debug','-echo','-fix-acf','-no-fsv3-st2fir','-no-fsv3-whiten'}, % ignore
        s.debug = 1;

      otherwise
        fprintf(2,'ERROR: Flag %s unrecognized\n',flag);
        s = [];
        return;

    end % --- switch(flag) ----- %

  end % while(narg <= ninputargs)

return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Check that there is at least one more argument %%
function arg1check(flag,nflag,nmax)
  if(nflag>nmax) 
    fprintf(1,'ERROR: Flag %s needs one argument',flag);
    error;
  end
return;
%--------------------------------------------------%
%% Check that there are at least two more arguments %%
function arg2check(flag,nflag,nmax)
  if(nflag > nmax-1 ) 
    fprintf(1,'ERROR: Flag %s needs two arguments',flag);
    error;
  end
return;


%--------------------------------------------------%
%% Check argument consistency, etc %%%
function s = check_params(s)

  fprintf(1,'Checking Parameters\n');

  s.nruns = size(s.invollist,1);
  npars = size(s.parlist,1);
  ntpxs = size(s.tpxlist,1);

  if(s.nruns < 1) 
    fprintf(2,'ERROR: No input volumes specified\n');
    s=[]; return;
  end

  if(s.nslices < 0)
    instem = deblank(s.invollist(1,:));
    if(~s.UseMRIread)
      [s.nslices nrows ncols ntrs] = fmri_bvoldim(instem);
    else
      mri = MRIread(instem,1);
      if(~isempty(mri)) s.nslices = mri.volsize(3);
      else              s.nslices = 0;
      end
    end
    if(s.nslices == 0) 
      fprintf(2,'ERROR: Volume %s does not exist\n',instem);
      s=[]; return;
    end      
    fprintf('nslices = %d\n',s.nslices);
  end

  if(npars ~= 0 & ~isempty(s.parname) ) 
    fprintf(2,'ERROR: Cannot specify both -p and -parname\n');
    s=[]; return;
  end

  if(npars == 0 & isempty(s.parname) ) 
    fprintf(2,'ERROR: No paradigm specified\n');
    s=[]; return;
  end

  if( ~isempty(s.parname) ) 
    for n = 1:s.nruns
      involpath = fast_dirname(deblank(s.invollist(n,:)));
      par = sprintf('%s/%s',involpath,s.parname);
      s.parlist = strvcat(s.parlist,par);
    end
    npars = size(s.parlist,1);
  end

  if(npars ~= s.nruns)
    fprintf(2,'ERROR: Number of input volumes (%d) and paradigms (%d) differ\n',...
                  s.nruns,npars);
    s=[]; return;
  end

  if(ntpxs ~= 0 & ntpxs ~= s.nruns)
    fprintf(2,'ERROR: Number of input volumes (%d) and tpexcl files (%d) differ\n',...
                  s.nruns,ntpxs);
    s=[]; return;
  end

  if(~isempty(s.extregfile) ) 
    for n = 1:s.nruns
      involpath = fast_dirname(deblank(s.invollist(n,:)));
      extregtmp = sprintf('%s/%s',involpath,s.extregfile);
      s.extreglist = strvcat(s.extreglist,extregtmp);
    end
  end

  if(size(s.hvol,1) ~= 1)
    fprintf(2,'ERROR: No output volume specified\n');
    s = []; return;
  end

  if(s.NoAutoWhiten) s.AutoWhiten = 0; end
  if(s.AutoWhiten)
    s.hvol = sprintf('%s0',s.hvol);
    fprintf('INFO: chaning output volume stem to %s for first stage\n');
  end
  if(s.AutoWhiten & s.WhitenFlag)
    fprintf('ERROR: cannot specify -autowhiten and -whiten\n');
    s = []; return;
  end
  if(s.WhitenFlag)
    if(size(s.WhtnMtxFile,1) ~= s.nruns)
      fprintf('ERROR: must spec nruns whtmtx files with -whiten\n');
      s = []; return;
    end
  end

  if(~isempty(s.HPF))
    if(s.HPF(1) < 0 | s.HPF(1) > 1 | s.HPF(2) < 0 | s.HPF(2) >= 1)
      fprintf(2,'ERROR: HPF Parameters out of range\n');
      s = []; return;
    end
  end

  if(length(s.TR) == 0)
    fprintf(2,'ERROR: No TR specified\n');
    s = []; return;
  end

  if(length(s.TotWin) == 0)
    fprintf(2,'ERROR: No Time Window specified \n');
    s = []; return;
    %fprintf(2,'INFO: No Time Window specified ...\n');
    %fprintf(2,' Setting to 20 sec\n');
    %s.TotWin = 20;
  end

  if(length(s.AcqOrder) > 0)
    if(~strcmpi(s.AcqOrder,'Linear') & ~strcmpi(s.AcqOrder,'Interleaved'))
     fprintf(2,'ERROR: Acquisition Order %s unknown (Linear or Interleaved)\n',...
              s.AcqOrder);
      s = []; return;
    end
  end

  if(length(s.TER) == 0) s.TER = s.TR; end

  if(s.firstslice < 0) 
    fprintf('ERROR: firstslice (%d) < 0',s.firstslice);
    s = []; return;
  end

  s.GammaFit = length(s.gfDelta);

  if(s.SaveErrCovMtx)
    s.ErrCovMtxStem = sprintf('%s-ecvm',s.hvol);
  end

  if(s.FWHM > 0 & s.HanRad > 0)
    fprintf('ERROR: Cannot specify both -hanrad and -fwhm\n');
    s = []; return;
  end

  if(s.FWHM > 0 & isempty(s.InPlaneRes))
    fprintf('ERROR: Need -ipr with -fwhm\n');
    s = []; return;
  end

  if(s.FWHM > 0 )
    s.HanRad = pi*s.FWHM/(2*s.InPlaneRes*acos(.5));
  end

  fprintf('AutoStimDur: %d\n',s.autostimdur);
  fprintf('StimDur: %g\n',s.stimdur);
  
  if(s.autostimdur & ~isempty(s.stimdur))
    fprintf('ERROR: cannot specify both autostimdur and stimdur\n');
    s = []; return;
  end
    
return;

%--------------------------------------------------%
%% Print data structure
function s = sxa_print_struct(s,fid)
  if(nargin == 1) fid = 1; end

  fprintf(fid,'Number of Runs: %d\n',s.nruns);

  fprintf(fid,'Input Volume List\n');
  for n = 1:size(s.invollist,1),
    fprintf(fid,'  %d  %s\n',n,s.invollist(n,:));    
  end

  fprintf(fid,'Input Pardigm File List\n');
  for n = 1:size(s.parlist,1),
    fprintf(fid,'  %d  %s\n',n,s.parlist(n,:));    
  end

  if(~isempty(s.tpxlist))
    fprintf(fid,'TP Exclude File List\n');
    for n = 1:size(s.tpxlist,1),
      fprintf(fid,'  %d  %s\n',n,s.tpxlist(n,:));    
    end
  end

  fprintf(fid,'Output Volume  %s\n',s.hvol);
  if(~isempty(s.betavol))
     fprintf(fid,'Beta Volume  %s\n',s.betavol);
  end
  if(~isempty(s.fomnibusvol))
    fprintf(fid,'F Omnibus Volume  %s\n',s.fomnibusvol);
  end
  if(~isempty(s.pomnibusvol))
    fprintf(fid,'Sig Omnibus Volume  %s\n',s.pomnibusvol);
  end
  fprintf(fid,'TR    %f\n',s.TR);
  fprintf(fid,'TER   %f\n',s.TER);
  fprintf(fid,'Total   Window  %g\n',s.TotWin);
  fprintf(fid,'PreStim Window  %g\n',s.PreStimWin);
  fprintf(fid,'Remove Baseline %d\n',s.RmBaseline);
  fprintf(fid,'Remove Trend    %d\n',s.RmTrend);
  fprintf(fid,'Remove QTrend   %d\n',s.QTrendFit);
  fprintf(fid,'Rescale Target  %g\n',s.RescaleTarget);
  fprintf(fid,'nSkip           %d\n',s.nSkip);
  fprintf(fid,'InPlane Res     %g\n',s.InPlaneRes);
  fprintf(fid,'FWHM            %g\n',s.FWHM);
  fprintf(fid,'Hanning Radius  %g\n',s.HanRad);
  fprintf(fid,'Time Offset     %g\n',s.TimeOffset);
  if(~isempty(s.AcqOrder))
    fprintf(fid,'Acquistion Order %s\n',s.AcqOrder);
  end

  fprintf(fid,'GammaFit        %d\n',s.GammaFit);
  for n = 1:s.GammaFit
    fprintf(fid,'%d  %g  %g\n',n,s.gfDelta,s.gfTau);
  end
  fprintf('GammaFit Alpha: %g\n',s.gfAlpha);
  fprintf('SPM HRF: %g\n',s.spmhrf);

  fprintf(fid,'Seg Brain/Air   %d\n',s.SegBrainAir);
  fprintf(fid,'SynthSeed       %d\n',s.SynthSeed);

  if(~isempty(s.ErrCovMtxStem))
    fprintf(fid,'ErrCovMtx Stem   %s\n',s.ErrCovMtxStem);
  end

  if(~isempty(s.WhtnMtxFile))
    fprintf(fid,'WhtnMtx File   %s\n',s.WhtnMtxFile);
  end

  if(~isempty(s.extregfile))
    fprintf(fid,'ExtReg File   %s\n',s.extregfile);
    fprintf(fid,'NExtReg       %d\n',s.nextreg);
    fprintf(fid,'ExtRegOrthog  %d\n',s.extregorthog);
  end

  fprintf(fid,'firstslice   %d\n',s.firstslice);
  fprintf(fid,'nslices      %d\n',s.nslices);
  fprintf(fid,'nyqreg       %d\n',s.nyqreg);

return;
%--------------------------------------------------%
function runlist = getrunlist(invollist)
  nruns = size(invollist,1);
  runlist = [];
  for run = 1:nruns
    invol = deblank(invollist(run,:));
    tmp = fast_dirname(invol);
    runid = basename(tmp);
    runno = sscanf(runid,'%d');
    runlist = [runlist runno];
  end
return;

