function r = fast_svdfunctcvm(varargin)
% r = fast_svdfunctcvm(varargin)


%
% fast_svdfunctcvm.m
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

version = 'fast_svdfunctcvm.m @FS_VERSION@';
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
s = check_params(s);
if(isempty(s)) return; end

sxa_print_struct(s,1);

TR  = s.TR;
TER = s.TER;
TW  = s.TotWin;
TPS = s.PreStimWin;
RmBaseline = s.RmBaseline;
RmTrend    = s.RmTrend;
GammaFit = s.GammaFit;
gfDelta = s.gfDelta;
gfTau =  s.gfTau;
nskip = s.nSkip;
TimeOffset = s.TimeOffset;
AcqOrder = s.AcqOrder;
parfilelist = s.parlist;
instemlist = s.invollist;
tpxlist = s.tpxlist;
svddirlist = s.svddirlist;
nruns = size(parfilelist,1);

%-------------------------------------------------%
Nfir = round(TW/TER);

% Load Whitening matrix %
if(~isempty(s.WhtnMtxFile))
   WhtnMtx = fmri_ldbfile(s.WhtnMtxFile);
   [wu ws wv] = svd(WhtnMtx);
   WhtnMtxSqrt = wu * sqrt(ws) * wv'; %'Whitening filter
else
   WhtnMtx = [];
end

%-- Determine Number of Conditions across all runs----------%
% Get list of unqiue condition IDs for all runs %
condlist = [];
for run = 1:nruns
  par = fmri_ldpar(deblank(parfilelist(run,:)));
  condid = par(:,2);
  clrun = unique(par(:,2));
  condlist = unique([clrun; condlist]);
end

% Remove -1 and 0 %
ind = find(condlist ~= -1 & condlist ~= 0);
condlist = condlist(ind);

Nnnc = length(condlist); % excludes null %
Nc = Nnnc + 1;

fprintf(1,'Conditions Found (%d): ',Nnnc);
fprintf(1,'%2d ',condlist);
fprintf(1,'\n');

% Check for holes in the list of condition numbers %
if(~isempty(find(diff(condlist)~=1)))
  fprintf(2,'fast_selxavg2: missing conditions\n');
  return;
end

% Count the number per condition %
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
end

SubSampRate = TR/TER;

% Get basic info from the first run %
instem = deblank(instemlist(1,:));
[nslices nrows ncols ntrs] = fmri_bvoldim(instem);

% Load the temporal covariance matrix %
tcvmfile = sprintf('%s.bfloat',s.tcvmstem);
tcvmall = fmri_ldbfile(tcvmfile);
if(isempty(tcvmall))
  fprintf('ERROR: loading %s\n',tcvmfile);
  return;
end
Myall = tcvmall;

%----------- Loop over jackknifed runs -----------------%
tic;
for jkrun = 1:nruns
  fprintf('Run = %d, %g -------------------------- \n',jkrun,toc);

  jkrunlist = find([1:nruns] ~= jkrun);

  %-----------------------------------------%
  % Build the desing matrix across all runs %
  X = [];
  nthrun = 1;
  for run = jkrunlist

    instem = deblank(instemlist(run,:));

    % Get number of TRs in this run %
    [nslices nrows ncols ntrs] = fmri_bvoldim(instem);
    nvoxs = nrows*ncols;

    % Time Point Exclusion %
    if(~isempty(tpxlist))
       TPExcludeFile = deblank(tpxlist(run,:));
    else
       TPExcludeFile = [];
    end
    [indTPExcl indTPIncl] = fast_ldtpexcl(TPExcludeFile,TR,ntrs,nskip);
    ntpx = length(indTPExcl);
    fprintf(1,'       Excluding %d Points: ',ntpx);
    fprintf(1,'%d ',indTPExcl);
    fprintf(1,'\n');

    % Create Baseline/Trend Components of Convolution Matrix %
    Xbaseline = []; Xtrend    = [];
    if(RmBaseline) Xbaseline = fast_baselinemtx(nthrun,ntrs,nruns-1); end
    if(RmTrend)    Xtrend    = fast_trendmtx(nthrun,ntrs,nruns-1); end

    % Load paradigm for this run %
    par = fmri_ldpar(deblank(parfilelist(run,:)));

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

    % Tranform for Fitting to Gamma Function(s) %
    if(GammaFit > 0)
      Xpar = fmri_scm2gcm(Xfir,Nnnc,TER,TPS,gfDelta,gfTau);
      Navgs_per_cond = GammaFit;
    else
      Xpar = Xfir;
      Navgs_per_cond = Nfir;
    end

    % Number of averages excluding the offsets and trends
    NTaskAvgs = Nnnc*Navgs_per_cond;

    % Create final Convolution Matrix %
    Xi = [Xpar Xbaseline Xtrend ];

    % Exlude Points %
    Xitmp = Xi;
    Xi(indTPExcl,:) = 0;

    % Whitening Filter %
    if(isempty(s.WhtnMtxFile)) 
      WhtnMtx = eye(ntrs); 
      WhtnMtxSqrt = eye(ntrs); 
    else
      [wr wc] = size(WhtnMtx);
      if(wr ~= ntrs)
        fprintf(2,'ERROR: Whitening Matrix is %d x %d, ntrs = %d\n',...
                wr,wc,ntrs);
        return;
      end
    end

    % Exlude Points %
    Xi(indTPExcl,:) = 0;

    X = [X; Xi];

    nthrun = nthrun + 1;
  end 
  % ------ Loop over runs to build design matrix -------- %

  My = fast_tcvm_rmrun(jkrun,Myall,nruns); % Jackknife
  T = X*inv(X'*X)*X'; 
  E = eye(size(T)) - T;
  Ms = T*My*T'; %'
  Me = E*My*E'; %'

  fname = sprintf('%s/x.bfloat',deblank(s.svddirlist(jkrun,:)));
  fmri_svbfile(X,fname);

  fname = sprintf('%s/t.bfloat',deblank(s.svddirlist(jkrun,:)));
  fmri_svbfile(T,fname);

  fname = sprintf('%s/e.bfloat',deblank(s.svddirlist(jkrun,:)));
  fmri_svbfile(E,fname);

  fprintf('Raw SVD \n',run);
  [Uy Sy dummy] = svd(My);
  fname = sprintf('%s/my.bfloat',deblank(s.svddirlist(jkrun,:)));
  fmri_svbfile(My,fname);
  fname = sprintf('%s/uy.bfloat',deblank(s.svddirlist(jkrun,:)));
  fmri_svbfile(Uy,fname);
  fname = sprintf('%s/sy.bfloat',deblank(s.svddirlist(jkrun,:)));
  fmri_svbfile(Sy,fname);

  fprintf('Signal SVD \n',run);
  [Us Ss dummy] = svd(Ms);
  fname = sprintf('%s/ms.bfloat',deblank(s.svddirlist(jkrun,:)));
  fmri_svbfile(Ms,fname);
  fname = sprintf('%s/us.bfloat',deblank(s.svddirlist(jkrun,:)));
  fmri_svbfile(Us,fname);
  fname = sprintf('%s/ss.bfloat',deblank(s.svddirlist(jkrun,:)));
  fmri_svbfile(Ss,fname);

  fprintf('ResErr SVD \n',run);
  [Ue Se dummy] = svd(Me);
  fname = sprintf('%s/me.bfloat',deblank(s.svddirlist(jkrun,:)));
  fmri_svbfile(Me,fname);
  fname = sprintf('%s/ue.bfloat',deblank(s.svddirlist(jkrun,:)));
  fmri_svbfile(Ue,fname);
  fname = sprintf('%s/se.bfloat',deblank(s.svddirlist(jkrun,:)));
  fmri_svbfile(Se,fname);

end

fprintf(1,'Done %g\n',toc);

r = 0;

return;
%---\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%-----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%


%------------- Print Usage ---------------------%
function print_usage(dummy)

  fprintf(1,'USAGE:\n');
  fprintf(1,'  fast_svdfunctcvm\n');
  fprintf(1,'     -i   invol ... \n');
  fprintf(1,'     -p   parfile ... \n');
  fprintf(1,'     -parname  parfile \n');
  fprintf(1,'     -tpx tpxfile ... \n');
  fprintf(1,'     -whtmtx whitening matrix file \n');
  fprintf(1,'     -svddir dir <-svddir dir>...\n');
  fprintf(1,'     -svdsubdir subdir \n');
  fprintf(1,'     -tcvm      stem \n');
  fprintf(1,'     -TR   TR\n');
  fprintf(1,'     -TER  TER\n');
  fprintf(1,'     -timewindow totwin  \n');
  fprintf(1,'     -prewindow  prewin  \n');
  fprintf(1,'     -nobaseline  \n');
  fprintf(1,'     -detrend  \n');
  fprintf(1,'     -nskip  n \n');
  fprintf(1,'     -gammafit delta tau \n');
  fprintf(1,'     -timeoffset t \n');
  fprintf(1,'     -acqorder  <linear or interleaved> \n');
  fprintf(1,'     -cfg        file \n');

return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = sxa_struct
  s.invollist      = '';
  s.parlist        = '';
  s.nruns          = 0;
  s.parname        = '';
  s.tpxlist        = '';
  s.WhtnMtxFile     = '';
  s.AutoWhiten     = 0;
  s.NoAutoWhiten   = 0;
  s.MinCond        = 1;
  s.hvol           = '';
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
  s.RescaleTarget = 0;
  s.nSkip  = 0;
  s.HanRad = 0;
  s.gfDelta = [];
  s.gfTau = [];
  s.TimeOffset = 0;
  s.AcqOrder = '';
  s.SynthSeed = 0;
  s.cfgfile = '';
  s.verbose = 0;
  s.firstslice = 0;
  s.nslices    = -1;
  s.eresdir    = '';
  s.sigestdir  = '';

  s.svdsubdir = '';
  s.svddirlist = '';
  s.tcvmstem   = '';

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

      case '-svddir',
        arg1check(flag,narg,ninputargs);
        s.svddirlist = strvcat(s.svddirlist,inputargs{narg});
        narg = narg + 1;

      case {'-svdsubdir'},
        arg1check(flag,narg,ninputargs);
        s.svdsubdir = inputargs{narg};
        narg = narg + 1;

      case {'-tcvm'},
        arg1check(flag,narg,ninputargs);
        s.tcvmstem = inputargs{narg};
        narg = narg + 1;

      case '-p',
        arg1check(flag,narg,ninputargs);
        s.parlist = strvcat(s.parlist,inputargs{narg});
        narg = narg + 1;

      case {'-parname'},
        arg1check(flag,narg,ninputargs);
        s.parname = inputargs{narg};
        narg = narg + 1;

      case {'-tpx','-tpexclfile'}
        arg1check(flag,narg,ninputargs);
        s.tpxlist = strvcat(s.tpxlist,inputargs{narg});
        narg = narg + 1;

      case {'-whtnmtx'}
        arg1check(flag,narg,ninputargs);
        s.WhtnMtxFile = inputargs{narg};
        narg = narg + 1;

      case {'-autowhiten'} % Arg is minimum condition for regularization
        arg1check(flag,narg,ninputargs);
        s.AutoWhiten = 1;
        s.MinCond = sscanf(inputargs{narg},'%f',1);
        s.SaveErrCovMtx = 1;
        narg = narg + 1;

      case {'-noautowhiten'} % To ease recursive calls
        s.NoAutoWhiten = 1;
        s.AutoWhiten = 0;
        s.MinCond = 1;

      case {'-o','-h'},
        arg1check(flag,narg,ninputargs);
        s.hvol = inputargs{narg};
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

      case {'-svecovmtx','-sverrcovmtx','-svecvm'},
        s.SaveErrCovMtx = 1;

      case {'-psc','percent'},
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
  
      case {'-rescale'}
        arg1check(flag,narg,ninputargs);
        s.RescaleTarget = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-nskip'}
        arg1check(flag,narg,ninputargs);
        s.nSkip = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-hanrad'}
        arg1check(flag,narg,ninputargs);
        s.HanRad = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-gammafit'}
        arg2check(flag,narg,ninputargs);
        gfDelta = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;
        gfTau   = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;
        s.gfDelta = [s.gfDelta gfDelta];
        s.gfTau   = [s.gfTau   gfTau];

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

      case '-verbose',
        s.verbose = 1;

      % ignore these guys %
      case {'-monly', '-nullcondid','-umask','-sveres','-svsignal'},
        arg1check(flag,narg,ninputargs);
        narg = narg + 1;

      case {'-debug','-echo'}, % ignore

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

  s.nruns  = size(s.invollist,1);
  npars    = size(s.parlist,1);
  nsvddirs = size(s.svddirlist,1);
  ntpxs    = size(s.tpxlist,1);

  if(isempty(s.tcvmstem)) 
    fprintf(2,'ERROR: No tcvm specified\n');
    s=[]; return;
  end

  if(s.nruns < 1) 
    fprintf(2,'ERROR: No input volumes specified\n');
    s=[]; return;
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


  if(nsvddirs ~= 0 & ~isempty(s.svdsubdir) ) 
    fprintf(2,'ERROR: Cannot specify both -svddir and -svdsubdir\n');
    s=[]; return;
  end

  if(nsvddirs == 0 & isempty(s.svdsubdir) ) 
    fprintf(2,'ERROR: No svddir specified\n');
    s=[]; return;
  end

  if( ~isempty(s.svdsubdir) ) 
    for n = 1:s.nruns
      involpath = fast_dirname(deblank(s.invollist(n,:)));
      svddir = sprintf('%s/%s',involpath,s.svdsubdir);
      s.svddirlist = strvcat(s.svddirlist,svddir);
    end
    nsvddirs = size(s.svddirlist,1);
  end

  if(nsvddirs ~= s.nruns)
    fprintf(2,'ERROR: Number of input volumes (%d) and svddirs (%d) differ\n',...
                  s.nruns,nsvddirs);
    s=[]; return;
  end

  if(ntpxs ~= 0 & ntpxs ~= s.nruns)
    fprintf(2,'ERROR: Number of input volumes (%d) and tpexcl files (%d) differ\n',...
                  s.nruns,ntpxs);
    s=[]; return;
  end

  if(s.NoAutoWhiten)
    s.AutoWhiten = 0;
  end
  if(s.AutoWhiten)
    s.hvol = sprintf('%s0',s.hvol);
    fprintf('INFO: chaning output volume stem to %s for first stage\n');
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

  s.GammaFit = length(s.gfDelta);

  if(s.SaveErrCovMtx)
    s.ErrCovMtxStem = sprintf('%s-ecvm',s.hvol);
  end

return;

%--------------------------------------------------%
%% Print data structure
function s = sxa_print_struct(s,fid)
  if(nargin == 1) fid = 1; end

  fprintf(fid,'Number of Runs: %d\n',s.nruns);
  fprintf(fid,'tcvm %s\n',s.tcvmstem);

  fprintf(fid,'Input Volume List\n');
  for n = 1:size(s.invollist,1),
    fprintf(fid,'  %d  %s\n',n,s.invollist(n,:));    
  end

  fprintf(fid,'Input Pardigm File List\n');
  for n = 1:size(s.parlist,1),
    fprintf(fid,'  %d  %s\n',n,s.parlist(n,:));    
  end

  fprintf(fid,'SVD Dir List\n');
  for n = 1:size(s.svddirlist,1),
    fprintf(fid,'  %d  %s\n',n,s.svddirlist(n,:));    
  end

  if(~isempty(s.tpxlist))
    fprintf(fid,'TP Exclude File List\n');
    for n = 1:size(s.tpxlist,1),
      fprintf(fid,'  %d  %s\n',n,s.tpxlist(n,:));    
    end
  end

  fprintf(fid,'TR    %f\n',s.TR);
  fprintf(fid,'TER   %f\n',s.TER);
  fprintf(fid,'Total   Window  %g\n',s.TotWin);
  fprintf(fid,'PreStim Window  %g\n',s.PreStimWin);
  fprintf(fid,'Remove Baseline %d\n',s.RmBaseline);
  fprintf(fid,'Remove Trend    %d\n',s.RmTrend);
  fprintf(fid,'nSkip           %d\n',s.nSkip);
  fprintf(fid,'Time Offset     %g\n',s.TimeOffset);
  if(~isempty(s.AcqOrder))
    fprintf(fid,'Acquistion Order %s\n',s.AcqOrder);
  end

  fprintf(fid,'GammaFit        %d\n',s.GammaFit);
  for n = 1:s.GammaFit
    fprintf(fid,'%d  %g  %g\n',n,s.gfDelta,s.gfTau);
  end

  if(~isempty(s.WhtnMtxFile))
    fprintf(fid,'WhtnMtx File   %s\n',s.WhtnMtxFile);
  end

return;
%--------------------------------------------------%


