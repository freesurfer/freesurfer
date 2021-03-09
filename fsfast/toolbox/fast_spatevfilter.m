function r = fast_spatevfilter(varargin)
% r = fast_spatev(varargin)


%
% fast_spatevfilter.m
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

version = 'fast_spatevfilter.m @FS_VERSION@';
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

% Get basic info about the input volume %
[nslices nrows ncols ntrs] = fmri_bvoldim(s.involid);
nv = nslices*nrows*ncols;

if(s.nkeep < 0) s.nkeep = ntrs; end
indkeep = 1:s.nkeep;

tic;
fprintf('Loading Raw SEVals %g\n',toc);
fname = sprintf('%s.bfloat',s.rawsevalid);
Sy = fmri_ldbfile(fname);
if(isempty(Sy))
  fprintf(2,'ERROR reading %s\n',fname);
  return;
end

fprintf('Loading Signal SEVals %g\n',toc);
fname = sprintf('%s.bfloat',s.signalsevalid);
Ss = fmri_ldbfile(fname);
if(isempty(Ss))
  fprintf(2,'ERROR reading %s\n',fname);
  return;
end


%--------------------------------------------------------%
tic;
for slice = s.firstslice:s.lastslice
  fprintf('Slice %2d,  %g -------- \n',slice,toc);

  fprintf('Loading Raw Data %g\n',toc);
  fname = sprintf('%s_%03d.bshort',s.involid,slice);
  y = fmri_ldbfile(fname);
  if(isempty(y))
    fprintf(2,'ERROR reading %s\n',fname);
    return;
  end
  [nr nc nt] = size(y);
  y = reshape(y, [nr*nc nt])'; %'

  fprintf('Loading Raw SEVs %g\n',toc);
  fname = sprintf('%s_%03d.bfloat',s.rawsevectid,slice);
  Vy = fmri_ldbfile(fname);
  if(isempty(Vy))
    fprintf(2,'ERROR reading %s\n',fname);
    return;
  end
  [nr nc nt] = size(Vy);
  Vy = reshape(Vy, [nr*nc nt]);

  fprintf('Loading Signal SEVs %g\n',toc);
  fname = sprintf('%s_%03d.bfloat',s.signalsevectid,slice);
  Vs = fmri_ldbfile(fname);
  if(isempty(Vs))
    fprintf(2,'ERROR reading %s\n',fname);
    return;
  end
  [nr nc nt] = size(Vs);
  Vs = reshape(Vs, [nr*nc nt]);

  fprintf('Filtering %g\n',toc);
  yhat = ((y*Vy) * inv(Sy) * (Vy'*Vs) * Ss) *Vs'; 

  fprintf('Saving %g\n',toc);
  yhat = reshape(yhat', [nr nc nt]); %'
  fname = sprintf('%s_%03d.bshort',s.outvolid,slice);
  fmri_svbfile(yhat,fname);

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
  fprintf(1,'  fast_spatevfilter\n');
  fprintf(1,'     -o         stem \n');
  fprintf(1,'     -i         stem \n');
  fprintf(1,'     -rawsevect  stem \n');
  fprintf(1,'     -rawseval   stem \n');
  fprintf(1,'     -signalsevect stem \n');
  fprintf(1,'     -signalseval  stem \n');
  fprintf(1,'     -nkeep n \n');
  fprintf(1,'     -regmethod method \n');
  fprintf(1,'     -fs firstslice \n');
  fprintf(1,'     -ns nslices \n');
return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = spatevfilter_struct
  s.involid       = '';
  s.outvolid      = '';
  s.rawsevectid      = '';
  s.rawsevalid      = '';
  s.signalsevectid   = '';
  s.signalsevalid   = '';
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
  s = spatevfilter_struct;
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
        s.involid = inputargs{narg};
        narg = narg + 1;

      case '-o',
        arg1check(flag,narg,ninputargs);
        s.outvolid = inputargs{narg};
        narg = narg + 1;

      case '-rawsevect',
        arg1check(flag,narg,ninputargs);
        s.rawsevectid = inputargs{narg};
        narg = narg + 1;

      case '-rawseval',
        arg1check(flag,narg,ninputargs);
        s.rawsevalid = inputargs{narg};
        narg = narg + 1;

      case '-signalsevect',
        arg1check(flag,narg,ninputargs);
        s.signalsevectid = inputargs{narg};
        narg = narg + 1;

      case '-signalseval',
        arg1check(flag,narg,ninputargs);
        s.signalsevalid = inputargs{narg};
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

  if(isempty(s.involid))
    fprintf(2,'ERROR: No input specified\n');
    s=[]; return;
  end
   
  if(isempty(s.outvolid))
    fprintf(2,'ERROR: No output specified\n');
    s=[]; return;
  end
   
  if(isempty(s.rawsevectid))
    fprintf(2,'ERROR: No raw sev specified\n');
    s=[]; return;
  end
   
  if(isempty(s.signalsevectid))
    fprintf(2,'ERROR: No signal sev specified\n');
    s=[]; return;
  end

  if(s.nslices < 0)
    instem = deblank(s.involid);
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

  fprintf(fid,'Input  volume %s\n',s.involid);
  fprintf(fid,'Output volume %s\n',s.outvolid);
  fprintf(fid,'Raw SEVect    %s\n',s.rawsevectid);
  fprintf(fid,'Raw SEVal     %s\n',s.rawsevalid);
  fprintf(fid,'Signal SEVect %s\n',s.signalsevectid);
  fprintf(fid,'Signal SEVal  %s\n',s.signalsevectid);
  fprintf(fid,'Reg Method    %s\n',s.regmethod);
  fprintf(fid,'nkeep         %d\n',s.nkeep);
  fprintf(fid,'firstslice    %d\n',s.firstslice);
  fprintf(fid,'nslices       %d\n',s.nslices);

return;
%--------------------------------------------------%


