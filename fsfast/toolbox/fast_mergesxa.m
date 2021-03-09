function r = fast_mergesxa(varargin)
% r = fast_mergesxa(varargin)


%
% fast_mergesxa.m
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

version = 'fast_mergesxa.m @FS_VERSION@';
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

%% Load and check that the sxa.dat files are consistent %%
Nc = 0;
Ntp = 0;
NPC = [];
for n = 1:s.ninputs
  instem = deblank(s.insxalist(n,:));
  datfile = sprintf('%s.dat',instem);
  sxa = fmri_lddat3(datfile);
  if(isempty(sxa))
    fprintf('ERROR: cannot open %s\n',datfile);
    return;
  end

  if(n == 1)
    sxa0 = sxa;
    datfile0 = datfile;
    DOFtot = sxa.DOF;
    NPC = sxa.Npercond;
    Ntp = sxa.Ntp;
  else
    err = 0;
    if(sxa.Version ~= sxa0.Version)       err = 1; end
    if(sxa.TR ~= sxa0.TR)                 err = 2; end
    if(sxa.TER ~= sxa0.TER)               err = 3; end
    if(sxa.TimeWindow ~= sxa0.TimeWindow) err = 4; end
    if(sxa.TPreStim ~= sxa0.TPreStim)     err = 5; end
    if(sxa.Nh ~= sxa0.Nh)                 err = 6; end
    if(err)
      fprintf('ERROR (%d): %s is incompatible with %s\n',...
	      err,datfile0,datfile);
      return;
    end
    NPC(1) = NPC(1) + sxa.Npercond(1);
    NPC = [NPC; sxa.Npercond(2:length(sxa.Npercond))];
    Ntp = Ntp + sxa.Ntp;
  end

  Nc = Nc + sxa.Nnnc;
end

nAvgsTot = Nc * sxa0.Nh;
SumXtX  = zeros(nAvgsTot);
hCovMtx = zeros(nAvgsTot);

instem = deblank(s.insxalist(1,:));
[nslices nrows ncols nframes] = fmri_bvoldim(instem);
nvxs = nrows*ncols;

%-----------------------------------------------------------%
tic;
for slice = 0:nslices-1,
  fprintf('%2d ',slice);

  hall = [];
  evarsum = 0;
  evarall = [];
  hoffsum = 0;
  DOF = 0;
  hi1 = 1;
  for n = 1:s.ninputs
    instem = deblank(s.insxalist(n,:));

    infile = sprintf('%s_%03d.bfloat',instem,slice);
    [h evar sxa] = fast_ldsxabfile(infile);
    if(isempty(h))
      fprintf('ERROR: could not load %s\n',infile);
      return;
    end
    nframes = size(h,3);
    h = reshape(h, [nvxs nframes])'; %'
    evar = reshape(evar, [nvxs 1])'; %'    

    offsetfile = sprintf('%s-offset_%03d.bfloat',instem,slice);
    hoff = fmri_ldbfile(offsetfile);
    if(isempty(hoff))
      fprintf('ERROR: could not load %s\n',offsetfile);
      return;
    end

    hall = [hall; h];
    evarsum = evarsum + sxa.DOF * evar;
    hoffsum = hoffsum + sxa.DOF * hoff;
    DOF = DOF + sxa.DOF;
    evarall = [evarall; evar];

    hi2 = hi1 + sxa.Nnnc*sxa.Nh - 1;
    SumXtX(hi1:hi2,hi1:hi2) = sxa.SumXtX;
    hCovMtx(hi1:hi2,hi1:hi2) = sxa.hCovMtx;
    hi1 = hi2 + 1;

  end % Loop over inputs 

  % Update sxa struct (only needs to be done once, but ...) %
  sxaF         = sxa0;
  sxaF.DOF     = DOF;
  sxaF.SumXtX  = SumXtX;
  sxaF.hCovMtx = hCovMtx;
  sxaF.Nnnc    = Nc;
  sxaF.Nc      = Nc + 1;
  sxaF.Npercond = NPC;
  sxaF.Ntp      = Ntp;

  % Save sxa struct to dat file %
  fname = sprintf('%s.dat',s.outstem);
  fmri_svdat2(fname,sxaF);

  % Compute the final offset %
  hoffF = hoffsum/DOF;
  offsetfile = sprintf('%s-offset_%03d.bfloat',s.outstem,slice);
  fmri_svbfile(hoffF,offsetfile);

  % Compute the final residual error variances %
  evarF = evarsum/DOF;

  % Convert to selavg format %
  halltmp = [zeros(sxaF.Nh,nvxs); hall]; % Add zero for cond 0
  halltmp = reshape(halltmp,[sxaF.Nh (Nc+1) nvxs]);

  hstd = sqrt( (diag(hCovMtx).*diag(SumXtX)) * evarF);
  hstd = [repmat(sqrt(evarF), [sxaF.Nh 1]); hstd]; % Add estd for cond 0
  hstd = reshape(hstd,[sxaF.Nh (Nc+1) nvxs]);

  %--- Merge Averages and StdDevs ---%
  tmp = zeros(sxaF.Nh,2,Nc+1,nvxs);
  tmp(:,1,:,:) = halltmp;
  tmp(:,2,:,:) = hstd;
  tmp = reshape(tmp,[sxaF.Nh*2*(Nc+1) nrows ncols ]);
  tmp = permute(tmp,[2 3 1]);

  % Save in selxavg format %
  fname = sprintf('%s_%03d.bfloat',s.outstem,slice);
  %fprintf(1,'  Saving data to %s \n',fname);
  fmri_svbfile(tmp,fname);

end % Loop over slices
fprintf('\n');

fprintf(1,'Done %g\n',toc);


r = 0;

return;
%---\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%-----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%


%------------- Print Usage ---------------------%
function print_usage(dummy)

  fprintf(1,'USAGE:\n');
  fprintf(1,'  fast_mergesxa\n');
  fprintf(1,'     -i insxavol1 -i ... (or use -ilist)\n');
  fprintf(1,'     -nlist ninputs (needed for -ilist)\n');
  fprintf(1,'     -ilist insxavol1 insxavol2 ... \n');
  fprintf(1,'     -o outsxavol \n');

return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = main_struct
  s.insxalist      = '';
  s.nlist          = -1;
  s.ninputs        = 0;
  s.outstem       = '';
return;

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

      case '-i',
        arg1check(flag,narg,ninputargs);
        s.insxalist = strvcat(s.insxalist,inputargs{narg});
        narg = narg + 1;

      case '-nlist',
        arg1check(flag,narg,ninputargs);
        s.nlist = sscanf(inputargs{narg},'%d');
        narg = narg + 1;

      case '-ilist',
        if(s.nlist == 0)
          fprintf('ERROR: -nlist must be specified before -ilist\n');
          s = []; return;
        end
        for n = 1:s.nlist
          arg1check(flag,narg,ninputargs);
          s.insxalist = strvcat(s.insxalist,inputargs{narg});
          narg = narg + 1;
        end

      case '-o',
        arg1check(flag,narg,ninputargs);
        s.outstem = inputargs{narg};
        narg = narg + 1;

      case '-verbose',
        s.verbose = 1;

      % ignore these guys %
      case {'-monly', '-umask'},
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

  s.ninputs = size(s.insxalist,1);
  if(s.ninputs < 1) 
    fprintf(2,'ERROR: No input volumes specified\n');
    s=[]; return;
  end

  if(isempty(s.outstem))
    fprintf(2,'ERROR: No output volume specified\n');
    s=[]; return;
  end

return;

%--------------------------------------------------%
%% Print data structure
function s = main_print_struct(s,fid)
  if(nargin == 1) fid = 1; end

  fprintf(fid,'Number of Runs: %d\n',s.ninputs);
  fprintf(fid,'Input Volume List\n');
  for n = 1:s.ninputs,
    fprintf(fid,'  %d  %s\n',n,s.insxalist(n,:));    
  end
  fprintf(fid,'Output Volume: %s\n',s.outstem);

return;
%--------------------------------------------------%


