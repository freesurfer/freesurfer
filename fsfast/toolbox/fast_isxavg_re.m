function r = fast_isxavg_re(varargin)
% r = fast_isxavg_re(varargin)
% Random effects intersubect averaging


%
% fast_isxavg_re.m
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

ninvols = size(s.invols,1);
s.dof = ninvols - 1;
lastslice = s.firstslice + s.nslices - 1;

fprintf(1,'_______________________________________________\n');
fprintf(1,'Random Effect Averaging Parameters\n');
isxavg_re_print_struct(s,1);
fprintf(1,'_______________________________________________\n');

%% Contrast Matrix %%
if(~isempty(s.cmtxfile))
  %ContrastMtx_0 = load(s.cmtxfile);
  %ContrastMtx_0 = [];
  % for some reason, cannot compile load(s.cmtxfile); %
  tmpfile = s.cmtxfile;
  tmp = load(tmpfile);
  R = tmp.ContrastMtx_0;
else 
  R = [];
  % Make R=I later when dims are known
end

% Synthesize ? %
if(s.synth)
  Seed = sum(100*clock);
  randn('state',Seed);
  fprintf(1,'INFO: Synthsizing with zero-mean/unit-variance (%d)\n',Seed);
end

% ---------- Go through each slice -------------- %
firstpass = 1;
for slice = s.firstslice:lastslice
  fprintf(1,'%2d ',slice);

  % ------ Load each input volume -------------%
  for n = 1:ninvols
    invol = deblank(s.invols(n,:));

    %% Load the data %%
    fname = sprintf('%s_%03d.bfloat',invol,slice);
    datname = sprintf('%s.dat',invol);
    if(fast_fileexists(datname))
      % Load as a selxavg volume 
      [yi evar hdrdat] = fast_ldsxabfile(fname);
    else
      % Load as a simple volume 
      yi = fmri_ldbfile(fname);
      hdrdat = [];
    end

    % Check consistency with first volume %
    [nr nc nt] = size(yi);
    if(n == 1)
       nr0 = nr; nc0 = nc;  nt0 = nt;  hdrdat0 = hdrdat;
       nv = nr*nc;
       if(isempty(R)) R = eye(nt);
       else
         if(size(R,2) ~= nt)
           fprintf(1,'ERROR: %s is inconsistent with contrast mtx\n',invol);
           fprintf(1,' nt = %d, Cols of R = %d\n',nt,size(R,2));
           error;
         end
       end
       nRcols = size(R,1);
       z = zeros(nRcols,nv,ninvols);
    else
      if(nr ~= nr0 | nc ~= nc0 | nt ~= nt0 | ...
         (prod(size(hdrdat0)) ~= prod(size(hdrdat)) ))
        fprintf(1,'ERROR: %s is inconsistent with %s\n',...
              s.invols(1,:),invol);
	fprintf('\nIt may be that the analysis was redefined after\n');
	fprintf('one of the sessions had been processed. Make sure\n');
	fprintf('that each session has been processed using the\n');
	fprintf('same analysis and that the resampling has been\n');
	fprintf('done after running selxavg.\n\n');
        return;
      end
    end % if(n == 1)

    % Check consistency with contrast matrix %
    if(size(R,2) ~= nt)
      fprintf('\nERROR: contrast matrix size is inconsistent with\n');
      fprintf('the number of parameters in the analysis. Check\n');
      fprintf('whether the analysis was redifined without redefining\n');
      fprintf('the contrast\n\n');
      return;
    end

    % Percent Signal Change %
    if(s.pctsigch)
      % Load the offset slice %
      fname = sprintf('%s-offset_%03d.bfloat',invol,slice);
      hoffset = fmri_ldbfile(fname);
      iz = find(hoffset == 0);
      hoffset(iz) = 10^10;
      yi = yi ./ repmat(hoffset,[1 1 nt]);
    end

    % Reshape yi into array nt X nv %
    yi = reshape(yi, [nv nt])'; % '

    % Synthesize ? %
    if(s.synth) yi = randn(size(yi));  end

    % Compute contrast %
    zi = R*yi;

    % Truncate all values of the specified sign to 0
    if( ~isempty(s.trunc) )
      if( strcmpi(s.trunc,'pos') )
        ind = find(zi > 0);
      end
      if( strcmpi(s.trunc,'neg') )
        ind = find(zi < 0);
      end
      if( ~isempty(ind) ) zi(ind) = 0; end
    end

    % Add to an array of all inputs %
    z(:,:,n) = zi;

  end % for n = 1:ninvols
  %-------------- All the data is loaded in --------------------%

  %-- Average Data and compute std dev ---%
  if(s.jackknife)
    zavg = zeros(nRcols,nv);
    zstd = zeros(nRcols,nv);
    for m = 1:nRcols
      v = squeeze(z(m,:,:)); 
      if(nv > 1) v = v'; end %' % ninvols X nv 
      [vavg vstd] = fmri_jackknife(v);
      zavg(m,:) = vavg;
      zstd(m,:) = vstd;
    end
  else
    zavg = mean(z,3);
    zstd = std(z,[],3);
  end

  % Make sure std does not equal zero %
  indzstd0 = find(zstd == 0);
  zstd(indzstd0) = 10^(-10);

  % Compute stats and sigs %%
  t = zavg./(zstd/sqrt(ninvols)); % Not dof
  p = tTest(s.dof, t, 100);

  % put the zero back %
  zstd(indzstd0) = 0;

  indtlz = find(t<0);
  p(indtlz) = -p(indtlz) ;

  %% Save average volume %%
  if(~isempty(s.avgvol))
    tmp = zavg;
    tmp = reshape(tmp', [nr nc nRcols]); %'
    fname = sprintf('%s_%03d.bfloat',s.avgvol,slice);
    fmri_svbfile(tmp,fname);
  end

  %% Save stddev volume %%
  if(~isempty(s.stdvol))
    tmp = zstd;
    tmp = reshape(tmp', [nr nc nRcols]); %'
    fname = sprintf('%s_%03d.bfloat',s.stdvol,slice);
    fmri_svbfile(tmp,fname);
    fname = sprintf('%s.dof',s.stdvol);
    fid = fopen(fname,'w');
    fprintf(fid,'%d\n',s.dof);
    fclose(fid);
  end

  %% Save t volume %%
  if(~isempty(s.tvol))
    tmp = reshape(t', [nr nc nRcols]); %'
    fname = sprintf('%s_%03d.bfloat',s.tvol,slice);
    fmri_svbfile(tmp,fname);
  end

  %% Save sig volume %%
  if(~isempty(s.sigvol))
    % Convert to log10 %
    tmp = -log10(abs(p));
    tmp(indtlz) = -tmp(indtlz);
    tmp = reshape(tmp', [nr nc nRcols]); %'
    fname = sprintf('%s_%03d.bfloat',s.sigvol,slice);
    fmri_svbfile(tmp,fname);
  end

  %% Compute, Save minsig and iminsig  %%
  if(~isempty(s.minsigvol) | ~isempty(s.iminsigvol))

    %% Compute minsig and iminsig  %%
    [pmin ipmin] = min(abs(p),[],1);
    indpmin = sub2ind(size(p),ipmin,[1:nv]);
    nBonferroni = size(p,1);
    pmin = nBonferroni*p(indpmin);

    % Save minsig %
    if(~isempty(s.minsigvol))
      indplz = find(pmin<0);
      tmp = -log10(abs(pmin));
      tmp(indplz) = -tmp(indplz);
      tmp = reshape(tmp', [nr nc 1]); %'
      fname = sprintf('%s_%03d.bfloat',s.minsigvol,slice);
      fmri_svbfile(tmp,fname);
    end

    % Save iminsig %
    if(~isempty(s.iminsigvol))
      tmp = reshape(ipmin', [nr nc 1]); %'
      fname = sprintf('%s_%03d.bfloat',s.iminsigvol,slice);
      fmri_svbfile(tmp,fname);
    end

  end  %%%% if(~isempty(s.minsigvol) | ~isempty(s.iminsigvol)) %%%

  firstpass = 0;

  if(rem(slice,10) == 9) fprintf(1,'\n'); end

end %%% for slice = s.firstslice:lastslice %%%

if(~isempty(s.minsigvol) | ~isempty(s.iminsigvol))
  fprintf(1,'INFO: Bonferroni Correction: %d\n',nBonferroni);
end

%% Save average volume %%
if(~isempty(s.avgvol))
  fname = sprintf('%s.rfxdat',s.avgvol);
  fid = fopen(fname,'w');
  if(fid == -1) 
    msg = sprintf('Could not open %s for writing\n',fname);
    qoe(msg); error(msg);
  end
  fprintf(fid,'RandomEffectsAveraging\n');
  isxavg_re_print_struct(s,fid);
  fclose(fid);
end

r = 0;
fprintf(1,'fmri_isxavg_re: completed SUCCESSFULLY\n');

return;
%----------------------------------------------------------%
%----------------------------------------------------------%
%----------------------------------------------------------%
%----------------------------------------------------------%


%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = isxavg_re_struct;
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
        s.invols = strvcat(s.invols,inputargs{narg});
        narg = narg + 1;

      case {'-firstslice', '-fs'}
        arg1check(flag,narg,ninputargs);
        s.firstslice = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-nslices', '-ns'}
        arg1check(flag,narg,ninputargs);
        s.nslices = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-sig',
        arg1check(flag,narg,ninputargs);
        s.sigvol = inputargs{narg};
        narg = narg + 1;

      case '-trunc', 
        arg1check(flag,narg,ninputargs);
        s.trunc = inputargs{narg};
        narg = narg + 1;

      case '-minsig',
        arg1check(flag,narg,ninputargs);
        s.minsigvol = inputargs{narg};
        narg = narg + 1;

      case '-iminsig',
        arg1check(flag,narg,ninputargs);
        s.iminsigvol = inputargs{narg};
        narg = narg + 1;

      case '-t',
        arg1check(flag,narg,ninputargs);
        s.tvol = inputargs{narg};
        narg = narg + 1;

      case '-avg',
        arg1check(flag,narg,ninputargs);
        s.avgvol = inputargs{narg};
        narg = narg + 1;

      case '-std',
        arg1check(flag,narg,ninputargs);
        s.stdvol = inputargs{narg};
        narg = narg + 1;

      case '-cmtx',
        arg1check(flag,narg,ninputargs);
        s.cmtxfile = inputargs{narg};
        narg = narg + 1;

      case '-logfile',
        arg1check(flag,narg,ninputargs);
        s.logfile = inputargs{narg};
        narg = narg + 1;

      case '-monly', % ignore
        arg1check(flag,narg,ninputargs);
        narg = narg + 1;

      case '-jackknife',
        s.jackknife = 1;

      case '-pctsigch',
        s.pctsigch = 1;

      case '-nojackknife',
        s.jackknife = 0;

      case '-synth',
        s.synth = 1;

      case '-nolog',
        s.nolog = 1;
        s.logfid = 0;

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
  if(size(s.invols,1) < 2) 
    fprintf(1,'ERROR: must have at least 2 inputs\n');
    error;
  end

  ninvols = size(s.invols,1);
  for n = 1:ninvols
    invol = deblank(s.invols(n,:));

    s.involfmt = fast_getvolformat(invol);
    if(isempty(s.involfmt))
      fprintf(1,'ERROR: could not determine format of input volume ');
      fprintf(1,' %s \n',invol);
      fprintf(1,'Volume may not exist.\n');
      qoe;error;
    end

    switch(s.involfmt)
      case 'bfile'
        [nslices nrows ncols nt] = fmri_bvoldim(invol);
	if(nslices == 0) s=[];  return; end
        if(s.nslices == -1) s.nslices = nslices; end
      case 'minc'
        if(~fast_fileexists(invol))
          fprintf(1,'ERROR: input volume %s does not exist',invol);
          error;
        end
      otherwise
        fprintf(1,'ERROR: cannot handle format %s',s.involfmt);
        error;
    end % switch(s.involfmt)

  end % for n = 1:ninvols

  if(isempty(s.cmtxfile) )
    fprintf(1,'INFO: no contrast specified, assuming identity \n');
  end

  if(~isempty(s.cmtxfile))
    if(~fast_fileexists(s.cmtxfile))
      fprintf(1,'ERROR: %s does not exist \n',s.cmtxfile);
      error;
    end  
  end

  if(isempty(s.sigvol) & isempty(s.tvol) & isempty(s.avgvol))
    fprintf(1,'ERROR: no output specified\n');
    error;
  end

  if(~isempty(s.trunc))
    if( ~strcmpi(s.trunc,'pos') & ~strcmpi(s.trunc,'neg') )
      fprintf(1,'ERROR: truncation must be pos or negative\n');
      error;
    end
  end

return;

%--------------------------------------------------%
%% Print Usage 
function print_usage
  fprintf(1,'USAGE:\n');
  fprintf(1,'  fast_isxavg_re\n');
  fprintf(1,'     -i invol1 -i invol2 ... \n');
  fprintf(1,'     -cmtx contrast \n');
  fprintf(1,'     -trunc pos, neg \n');
  fprintf(1,'     -t       ttestvol \n');
  fprintf(1,'     -sig     sigvol \n');
  fprintf(1,'     -minsig  minsigvol \n');
  fprintf(1,'     -iminsig iminsigvol \n');
  fprintf(1,'     -avg     avgvol \n');
  fprintf(1,'     -std     stdvol \n');
  fprintf(1,'     -pctsigch \n');
  fprintf(1,'     -firstslice sliceno  \n');
  fprintf(1,'     -nslices    nslices  \n');
  fprintf(1,'     -nojackknife  \n');
return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = isxavg_re_struct
  s.invols    = '';
  s.involfmt  = '';
  s.involtype  = '';
  s.firstslice = 0;
  s.nslices    = -1;
  s.sigvol     = '';
  s.minsigvol  = '';
  s.iminsigvol = '';
  s.tvol   = '';
  s.avgvol = '';
  s.stdvol = '';
  s.pctsigch = 0;
  s.cmtxfile = '';
  s.cmtx     = []; 
  s.jackknife = 1;
  s.logfile = '';
  s.logfid  = 1;
  s.nolog   = 0;
  s.synth   = 0;
  s.dof     = 0;
  s.trunc   = '';
return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Print data structure
function s = isxavg_re_print_struct(s,fid)
  if(nargin == 1) fid = 1; end

  ninvols = size(s.invols,1);
  fprintf(fid,'ninvols    %d\n',ninvols);
  for n = 1:ninvols
    fprintf(fid,'  invol      %s\n',s.invols(n,:));
  end
  fprintf(fid,'involfmt    %s\n',s.involfmt);
  fprintf(fid,'involtype   %s\n',s.involtype);
  fprintf(fid,'firstslice  %d\n',s.firstslice);
  fprintf(fid,'nslices     %d\n',s.nslices);
  fprintf(fid,'sigvol     %s\n',s.sigvol);
  fprintf(fid,'tvol       %s\n',s.tvol);
  fprintf(fid,'avgvol     %s\n',s.avgvol);
  fprintf(fid,'stdvol     %s\n',s.stdvol);
  fprintf(fid,'pctsigch   %d\n',s.pctsigch);
  fprintf(fid,'cmtxfile    %s\n',s.cmtxfile);
  fprintf(fid,'jackknife   %d\n',s.jackknife);
  fprintf(fid,'logfile     %s\n',s.logfile);
  fprintf(fid,'nolog       %d\n',s.nolog);
  fprintf(fid,'dof         %d\n',s.dof);

return;
%--------------------------------------------------%

