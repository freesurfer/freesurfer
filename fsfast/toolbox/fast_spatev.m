function r = fast_spatev(varargin)
% r = fast_spatev(varargin)


%
% fast_spatev.m
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

version = 'fast_spatev.m @FS_VERSION@';
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

% Get basic info from the first run %
instem = deblank(s.invollist(1,:));
[nslices nrows ncols ntrs] = fmri_bvoldim(instem);
nv = nslices*nrows*ncols;

if(s.nkeep < 0) s.nkeep = ntrs; end
indkeep = 1:s.nkeep;

% Read in Jackknifed Design Matrix %
fname = sprintf('%s/x.bfloat',s.svddir);
X = fmri_ldbfile(fname);
if(isempty(X))
  fprintf(2,'ERROR reading %s\n',fname);
  return;
end
T = X*inv(X'*X)*X';
E = eye(size(T)) - T;

% Read in Jackknifed Raw EVects %
fname = sprintf('%s/uy.bfloat',s.svddir);
Uy = fmri_ldbfile(fname);
if(isempty(Uy))
  fprintf(2,'ERROR reading %s\n',fname);
  return;
end
% Read in Jackknifed Raw EVals %
fname = sprintf('%s/sy.bfloat',s.svddir);
Sy = fmri_ldbfile(fname);
if(isempty(Sy))
  fprintf(2,'ERROR reading %s\n',fname);
  return;
end

% Read in Jackknifed Signal EVects %
fname = sprintf('%s/us.bfloat',s.svddir);
Us = fmri_ldbfile(fname);
if(isempty(Us))
  fprintf(2,'ERROR reading %s\n',fname);
  return;
end
% Read in Jackknifed Signal EVals %
fname = sprintf('%s/ss.bfloat',s.svddir);
Ss = fmri_ldbfile(fname);
if(isempty(Ss))
  fprintf(2,'ERROR reading %s\n',fname);
  return;
end

% Read in Jackknifed ResErr EVects %
fname = sprintf('%s/ue.bfloat',s.svddir);
Ue = fmri_ldbfile(fname);
if(isempty(Ue))
  fprintf(2,'ERROR reading %s\n',fname);
  return;
end
% Read in Jackknifed ResErr EVals %
fname = sprintf('%s/se.bfloat',s.svddir);
Se = fmri_ldbfile(fname);
if(isempty(Se))
  fprintf(2,'ERROR reading %s\n',fname);
  return;
end

% Compute Qs %
Qy = Uy*inv(sqrt(Sy*nv));
Qs = Us*inv(sqrt(Ss*nv));
Qe = Ue*inv(sqrt(Se*nv));

Qy = Qy(:,indkeep);
Qs = Qs(:,indkeep);
Qe = Qe(:,indkeep);

%--------------------------------------------------------%
tic;
for slice = s.firstslice:s.lastslice
  fprintf('Slice %2d,  %g -------- \n',slice,toc);

  yall = [];
  for run = 1:s.nruns
    fprintf('Loading Run %d, %g\n',run,toc);
    instem = deblank(s.invollist(run,:));
    fname = sprintf('%s_%03d.bshort',instem,slice);
    y = fmri_ldbfile(fname);
    if(isempty(y))
      fprintf(2,'ERROR reading %s\n',fname);
      return;
    end
    [nr nc nt] = size(y);
    y = reshape(y, [nr*nc nt])'; %'
    yall = [yall; y];
  end

  fprintf('Computing Raw SEV %g\n',toc);
  Vy = yall' * Qy; %'
  Vy = reshape(Vy, [nrows ncols s.nkeep]);
  fname = sprintf('%s-vy_%03d.bfloat',s.sevstem,slice);
  fmri_svbfile(Vy,fname);

  fprintf('Computing Signal SEV %g\n',toc);
  Vs = (T*yall)' * Qs; %'
  % save sv Vs T yall Qs Us Ss;
  Vs = reshape(Vs, [nrows ncols s.nkeep]);
  fname = sprintf('%s-vs_%03d.bfloat',s.sevstem,slice);
  fmri_svbfile(Vs,fname);


  fprintf('Computing ResErr SEV %g\n',toc);
  Ve = (E*yall)' * Qe;  %'
  Ve = reshape(Ve, [nrows ncols s.nkeep]); 
  fname = sprintf('%s-ve_%03d.bfloat',s.sevstem,slice);
  fmri_svbfile(Ve,fname);

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
  fprintf(1,'  fast_spatev\n');
  fprintf(1,'     -svddir dir\n');
  fprintf(1,'     -i   invol ... \n');
  fprintf(1,'     -sev stem \n');
  fprintf(1,'     -nkeep n \n');
  fprintf(1,'     -regmethod method \n');
  fprintf(1,'     -fs firstslice \n');
  fprintf(1,'     -ns nslices \n');
return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = spatev_struct
  s.invollist      = '';
  s.nruns          = 0;
  s.svddir         = '';
  s.sevstem        = '';
  s.nkeep          = -1;
  s.regmethod      = '';
  s.verbose = 0;
  s.firstslice = 0;
  s.nslices    = -1;
return;

%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = spatev_struct;
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
        s.invollist = strvcat(s.invollist,inputargs{narg});
        narg = narg + 1;

      case {'-svddir'},
        arg1check(flag,narg,ninputargs);
        s.svddir = inputargs{narg};
        narg = narg + 1;

      case {'-sev'}
        arg1check(flag,narg,ninputargs);
        s.sevstem = inputargs{narg};
        narg = narg + 1;

      case {'-regmethod'}
        arg1check(flag,narg,ninputargs);
        s.regmethod = inputargs{narg};
        narg = narg + 1;

      case {'-nkeep'}
        arg1check(flag,narg,ninputargs);
        s.nkeep = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-fs','-firstslice'}
        arg1check(flag,narg,ninputargs);
        s.firstslice = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-ns','-nslices'}
        arg1check(flag,narg,ninputargs);
        s.nslices = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-verbose',
        s.verbose = 1;

      % ignore these guys %
      case {'-monly','umask'},
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
%% Check argument consistency, etc %%%
function s = check_params(s)

  fprintf(1,'Checking Parameters\n');

  s.nruns = size(s.invollist,1);

  if(isempty(s.svddir))
    fprintf(2,'ERROR: No svd dir specified\n');
    s=[]; return;
  end
   
  if(isempty(s.sevstem))
    fprintf(2,'ERROR: No sev stem specified\n');
    s=[]; return;
  end

  if(s.nruns < 1) 
    fprintf(2,'ERROR: No input volumes specified\n');
    s=[]; return;
  end

  if(s.nslices < 0)
    instem = deblank(s.invollist(1,:));
    [s.nslices nrows ncols ntrs] = fmri_bvoldim(instem);
  end

  if(s.firstslice < 0) 
    msg = sprintf('ERROR: firstslice (%d) < 0',s.firstslice);
    s = []; return;
  end

  s.lastslice = s.firstslice + s.nslices - 1;

return;

%--------------------------------------------------%
%% Print data structure
function s = sxa_print_struct(s,fid)
  if(nargin == 1) fid = 1; end

  fprintf(fid,'svd dir %s\n',s.svddir);
  fprintf(fid,'Number of Runs: %d\n',s.nruns);
  fprintf(fid,'Input Volume List\n');
  for n = 1:size(s.invollist,1),
    fprintf(fid,'  %d  %s\n',n,s.invollist(n,:));    
  end

  fprintf(fid,'Spat EV Volume  %s\n',s.sevstem);
  fprintf(fid,'Reg Method   %s\n',s.regmethod);
  fprintf(fid,'nkeep        %d\n',s.nkeep);
  fprintf(fid,'firstslice   %d\n',s.firstslice);
  fprintf(fid,'nslices      %d\n',s.nslices);

return;
%--------------------------------------------------%


