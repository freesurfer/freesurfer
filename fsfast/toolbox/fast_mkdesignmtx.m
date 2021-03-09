function r = fast_mkdesignmtx(varargin)
% r = fast_mkdesignmtx(varargin)


%
% fast_mkdesignmtx.m
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

version = 'fast_mkdesignmtx.m @FS_VERSION@';
fprintf(1,'%s\n',version);
r = 1;

%% Print usage if there are no arguments %%
if(nargin == 0)
  print_usage;
  return;
end

%% Parse the arguments %%
s = parse_args(varargin);
if(isempty(s)) return; end
s = check_params(s);
if(isempty(s)) return; end

if(s.verbose) print_main_struct(s,1); end


%-- Determine Number of Conditions across all runs----------%
% Get list of unqiue condition IDs for all runs %
condlist = [];
for run = 1:s.nruns
  par = fmri_ldpar(deblank(s.parlist(run,:)));
  condid = par(:,2);
  clrun = unique(par(:,2));
  condlist = unique([clrun; condlist]);
end

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

SubSampRate = s.TR/s.TER;
Nfir = round(s.TotWin/s.TER);

Xall = [];
for run = 1:s.nruns
  par = fmri_ldpar(deblank(s.parlist(run,:)));

  % Convert paradigm to FIR stimulus convolution matrix %
  Xfir = fmri_par2scm(par,Nc,SubSampRate*s.ntrs,s.TER,Nfir,s.PreStimWin);

  % For Sub-TR Estimation %
  if(s.TR ~= s.TER)
    Xfirtmp = Xfir;
    nn = [1:SubSampRate:size(Xfirtmp,1)];
    Xfir = Xfirtmp(nn,:);
  end

  % Tranform for Fitting to Gamma Function(s) %
  if(s.GammaFit > 0)
    Xpar = fmri_scm2gcm(Xfir,Nnnc,s.TER,s.PreStimWin,s.gfDelta,s.gfTau);
    Navgs_per_cond = s.GammaFit;
  else
    Xpar = Xfir;
    Navgs_per_cond = Nfir;
  end

  % Create Baseline/Trend Components of Convolution Matrix %
  Xmean = []; Xtrend    = [];
  if(s.MeanFit)  Xmean  = fast_baselinemtx(run,s.ntrs,s.nruns); end
  if(s.TrendFit) Xtrend = fast_trendmtx(run,s.ntrs,s.nruns); end

  % Create final Convolution Matrix %
  Xi = [Xpar Xmean Xtrend ];

  % Vertically (temporally) concatenate %
  Xall = [Xall; Xi];

  if(~isempty(s.xlist))
    X = Xpar; % Do not save detrending components
    fname = deblank(s.xlist(run,:));
    save(fname,'-v4','X');
  end

end

if(~isempty(s.xall)) save(s.xall,'-v4','Xall'); end

XtX = Xall'*Xall; %'
iXtX = inv(XtX);

%-----------------------------------------------------------%
% Compute the matrix which extracts the task-related signal %
Ntotavgs  = size(XtX,1);
Ntaskavgs = Navgs_per_cond * Nnnc;
nn = zeros(Ntotavgs,1);
nn(1:Ntaskavgs) = 1;
R = diag(nn);

%--------------------------------------------------------%
% Eall is a filter-type matrix which, when multiplied by 
% the raw input will result in the residual error. Note
% that it is (Ntrs*Nruns)X(Ntrs*Nruns) and cannot be 
% decomposed into runs.
Eall = eye(s.nruns*s.ntrs) - Xall * iXtX * Xall'; %'
if(~isempty(s.eall)) save(s.eall,'-v4','Eall'); end

%--------------------------------------------------------%
% Sall is a filter-type matrix which, when multiplied by 
% the raw input will result in the task-related signal. 
% Note that it is (Ntrs*Nruns)X(Ntrs*Nruns) and cannot be 
% decomposed into runs.
Sall =  Xall * R * iXtX * Xall'; %'
if(~isempty(s.sall)) save(s.sall,'-v4','Sall'); end

%-----------------------------------------------------------------%
% Compute variance reduction factor of the task-related components
iXtXTask = iXtX(1:Ntaskavgs,1:Ntaskavgs);
d = 1./diag(iXtXTask);
VRFMean = mean(d);
VRFMin  = min(d);
VRFMax  = max(d);
VRFStd  = std(d);

fprintf('\n');
fprintf('Variance Reduction Factor Stats:\n');
fprintf('VRF Mean: %g\n',VRFMean);
fprintf('VRF Min:  %g\n',VRFMin);
fprintf('VRF Max:  %g\n',VRFMax);
fprintf('VRF Std:  %g\n',VRFStd);
fprintf('\n');

r = 0;
return;
%---\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%-----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%


%------------- Print Usage ---------------------%
function print_usage(dummy)

  fprintf(1,'USAGE:\n');
  fprintf(1,'  fast_mkdesignmtx\n');
  fprintf(1,'     -p    parfile1 ... \n');
  fprintf(1,'     -ntrs n\n');
  fprintf(1,'     -TR         TR\n');
  fprintf(1,'     -TER        TER\n');
  fprintf(1,'     -timewindow totwin  \n');
  fprintf(1,'     -prewindow  prewin  \n');
  fprintf(1,'     -nomeanfit    \n');
  fprintf(1,'     -notrendfit   \n');
  fprintf(1,'     -gammafit delta tau \n');
  fprintf(1,'\n');
  fprintf(1,'     -x x1 <-x x2> ... \n');
  fprintf(1,'     -xall xallfile \n');
  fprintf(1,'     -eall eallfile \n');
  fprintf(1,'     -sall sallfile \n');
  fprintf(1,'\n');
return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = main_struct
  s.parlist        = '';
  s.nconditions    = 0;
  s.nruns          = 0;
  s.ntrs           = [];
  s.TR    = '';
  s.TER    = '';
  s.TotWin      = '';
  s.PreStimWin  = 0;
  s.MeanFit  = 1;
  s.TrendFit = 1;
  s.gfDelta = [];
  s.gfTau = [];

  s.xlist = [];
  s.xall  = [];
  s.eall = [];
  s.sall = [];

  s.verbose = 0;
return;

%--------------------------------------------------%
%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = main_struct;
  inputargs = varargin{1};
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

      case '-p',
        arg1check(flag,narg,ninputargs);
        s.parlist = strvcat(s.parlist,inputargs{narg});
        narg = narg + 1;

      case {'-ntrs','-ntp'}
        arg1check(flag,narg,ninputargs);
        s.ntrs = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-nconditions','-nc'}
        arg1check(flag,narg,ninputargs);
        s.nconditions = sscanf(inputargs{narg},'%d',1);
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

      case {'-prewindow','-prewin','-prestim','-tprestim'}
        arg1check(flag,narg,ninputargs);
        s.PreStimWin = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-meanfit','-fitmean','-baseline'}
        s.MeanFit = 1;
  
      case {'-nomeanfit','-nobaseline'}
        s.MeanFit = 0;
  
      case {'-trendfit','-fittrend','-detrend'}
        s.TrendFit = 1;

      case {'-notrendfit','-nodetrend'}
        s.TrendFit = 0;
  
      case {'-gammafit'}
        arg2check(flag,narg,ninputargs);
        gfDelta = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;
        gfTau   = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;
        s.gfDelta = [s.gfDelta gfDelta];
        s.gfTau   = [s.gfTau   gfTau];

      case '-x',
        arg1check(flag,narg,ninputargs);
        s.xlist = strvcat(s.xlist,inputargs{narg});
        narg = narg + 1;

      case '-xall',
        arg1check(flag,narg,ninputargs);
        s.xall = inputargs{narg};
        narg = narg + 1;

      case '-eall',
        arg1check(flag,narg,ninputargs);
        s.eall = inputargs{narg};
        narg = narg + 1;

      case '-sall',
        arg1check(flag,narg,ninputargs);
        s.sall = inputargs{narg};
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

  s.nruns = size(s.parlist,1);

  if(s.nruns < 1) 
    fprintf(2,'ERROR: No input paradigms specified\n');
    s=[]; return;
  end

  if(isempty(s.ntrs) )
    fprintf(2,'ERROR: must specify the number of TRs per run\n');
    s=[]; return;
  end

  nxlist = size(s.xlist,1);
  if(nxlist ~= 0 & nxlist  ~= s.nruns ) 
    fprintf(2,'ERROR: length of xlist (%d) does not equal number of runs (%d)\n',nxlist,nruns);
    s=[]; return;
  end

  if(length(s.TR) == 0)
    fprintf(2,'ERROR: No TR specified\n');
    s = []; return;
  end

  if(length(s.TotWin) == 0)
    fprintf(2,'ERROR: No Time Window specified \n');
    s = []; return;
  end

  if(length(s.TER) == 0) s.TER = s.TR; end

  s.GammaFit = length(s.gfDelta);

return;

%--------------------------------------------------%
%% Print data structure
function s = print_main_struct(s,fid)
  if(nargin == 1) fid = 1; end

  fprintf(fid,'Number of Runs: %d\n',s.nruns);

  fprintf(fid,'Input Pardigm File List\n');
  for n = 1:size(s.parlist,1),
    fprintf(fid,'  %d  %s\n',n,s.parlist(n,:));    
    if(~isempty(s.xlist)) fprintf(fid,'  %d  %s\n',n,s.xlist(n,:)); end
  end

  fprintf(fid,'nTRs  %d\n',s.ntrs);
  fprintf(fid,'TR    %f\n',s.TR);
  fprintf(fid,'TER   %f\n',s.TER);
  fprintf(fid,'Total   Window  %g\n',s.TotWin);
  fprintf(fid,'PreStim Window  %g\n',s.PreStimWin);
  fprintf(fid,'Remove Baseline %d\n',s.MeanFit);
  fprintf(fid,'Remove Trend    %d\n',s.TrendFit);
  fprintf(fid,'GammaFit        %d\n',s.GammaFit);
  for n = 1:s.GammaFit
    fprintf(fid,'%d  %g  %g\n',n,s.gfDelta,s.gfTau);
  end
return;
%--------------------------------------------------%


