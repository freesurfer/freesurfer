function r = fast_evspatfil(varargin)
% r = fast_evspatfil(varargin)
%
% Spatially filters by projecting out a given set
% of spatial eigenvectors as computed by fast_evfunc.m
%
%


%
% fast_evspatfil.m
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

version = 'fast_evspatfil.m @FS_VERSION@';
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
dump_params(s);

% Get basic info from the first run %
instem = deblank(s.involid(1,:));
[nslices nrows ncols ntrs] = fmri_bvoldim(instem);
nv = nslices*nrows*ncols;
nvslice = nrows*ncols;

[nslices nrows ncols ntrs endian bext] = fmri_bvoldim(s.involid);
if(nslices == 0)
  fprintf('ERROR: reading %s\n',s.involid);
  return;
end

% Set output extension to input extension %
if(isempty(s.obext)) s.obext = bext; end

tic;

%----------------- brain masking ------------------------------%
if(~isempty(s.maskid))
  s.mask = fmri_ldbvolume(s.maskid);
  if(isempty(s.mask)) 
    fprintf('ERROR: could not load %s\n',s.maskid);
    return;
  end 
  s.mask = s.mask(:,:,:,1);
  s.mask = s.mask > s.maskthresh;
elseif(~isempty(s.automaskthresh))
  fprintf('Computing AutoMask\n');
  instem = deblank(s.involid);
  f = fmri_ldbvolume(instem);
  f = f(:,:,:,1);
  fmn = mean(reshape1d(f));
  absthresh = fmn * s.automaskthresh;
  iover = find(f > absthresh);
  nover = length(iover);
  fprintf('  nover = %d, globalmn = %g, overmn = %g\n',...
          nover,fmn,mean(f(iover)));

  if(nover == 0)
    fprintf('ERROR: could not find any voxels above threshold\n');
    return;
  end 
  s.mask = zeros(size(f));
  s.mask(iover) = 1;

  if(~isempty(s.automaskid))
    fprintf('  Saving automask to %s\n',s.automaskid);
    fmri_svbvolume(s.mask,s.automaskid);
  end
else
  fprintf('INFO: not using a mask\n');
  s.mask = ones(nslices,nrows,ncols);
end 
if(s.cutends)
  fprintf('INFO: cutting off ends\n');
  s.mask(1,:,:) = 0;
  s.mask(nslices,:,:) = 0;
end

nmasktot = length(find(s.mask == 1));
if(nmasktot == 0)
  fprintf('ERROR: no voxels found in mask\n');
  return;
end 
fprintf('INFO: found %d voxels in mask\n',nmasktot);

%----------------- done setting up mask ----------------------------------%

%----------------- projection masking ------------------------------%
if(~isempty(s.projmaskid))
  fprintf('Loading Projection Mask %s\n',s.projmaskid);
  s.projmask = fmri_ldbvolume(s.projmaskid);
  if(isempty(s.projmask)) 
    fprintf('ERROR: could not load %s\n',s.projmaskid);
    return;
  end 
  fprintf('Projection Mask Threshold = %g\n',s.projmaskthresh);
  s.projmask = s.projmask(:,:,:,1);
  s.projmask = abs(s.projmask) > s.projmaskthresh;
  nprojmask = length(find(s.projmask==1));
  fprintf('Projection Mask has %d voxels\n',nprojmask);
else
  s.projmask = [];
end

% ------ Construct Drift Removal Matrix ----------- %
if(s.pforder > 0)
  X = fast_polytrendmtx(1,ntrs,1,s.pforder);
  D = eye(ntrs) - X*inv(X'*X)*X';
else
  D = [];
end

%--- Determine which components to keep -----------------------%
%indkeep = 1:s.nfilter;

%--------------------------------------------------------------%
sev_sumsqrd = 0;
fprintf('Starting Projection Loop %g\n',toc);
Pfil = 0;
Ppve = 0;
TCVMsum  = 0;
for sliceno = s.firstslice:s.lastslice
  fprintf('%3d ',sliceno);

  imask = find(s.mask(sliceno+1,:,:)==1);
  nmask = length(imask);
  if(nmask > 0)

    % Load the spat ev slice %
    sevstem = deblank(s.sevvolid);
    [nslices nrows ncols nevs endian bext] = fmri_bvoldim(sevstem);
    if(nslices == 0)
       fprintf('ERROR: could not load %s\n',sevstem);
       return;
    end
    fname = sprintf('%s_%03d.%s',sevstem,sliceno,bext);
    sev_all  = fmri_ldbfile(fname);
    sev_all  = reshape(sev_all,[nvslice nevs]); % Dont transpose

    % Exclude voxels inv sev from the projection mask %
    % This is NOT to be done on the filtering loop %
    if(~isempty(s.projmask))
      iprojmask = find(s.projmask(sliceno+1,:,:)==1);
      sev_all(iprojmask,:) = 0;
    end

    % Use the masked sev in this stage; the filter stage uses unmasked
    sev_all = sev_all(imask,:);

    % Zero out the sevs to be excluded
    if(~isempty(s.sevexcl)) sev_all(:,s.sevexcl) = 0; end

    % Compute Sum of the Squares for each vector. This is 
    % necessary to correct the PVE when a mask is used 
    sev_sumsqrd = sev_sumsqrd + sum(sev_all.^2)'; %'

    if(~isempty(s.nfilter)) sev = sev_all(:,1:s.nfilter); 
    else                    sev = sev_all;
    end

    % Load the input slice %
    instem = deblank(s.involid);
    [nslices nrows ncols ntrs endian bext] = fmri_bvoldim(instem);
    if(nslices == 0)
       fprintf('ERROR: could not load %s\n',instem);
       return;
    end
    fname = sprintf('%s_%03d.%s',instem,sliceno,bext);
    yrun = fmri_ldbfile(fname);
    yrun = reshape(yrun,[nvslice ntrs])'; %'
    yrun = yrun(:,imask);

    if(~isempty(s.outvolid))
      % Pfil is Ntp X Nfilter. Each column is an estimate of the 
      % noise temporal eigenvector
      Pfil = Pfil + yrun*sev;
    end

    if(~isempty(s.pvefile) | ~isempty(s.projtcid) )
      % Compute TCVM of the "raw" and projected data %
      if(~isempty(D)) yrun = D*yrun; end % Remove mean and/or trend 
      TCVMslice = yrun*yrun';%'
      TCVMsum = TCVMsum + TCVMslice;
      Ppve = Ppve + yrun*sev_all;
    end 

  end %if(nmask > 0)%

end % Loop over slices
fprintf('\n');

%--------------------------------------------------------------%
if(~isempty(s.pvefile))
  TCVM  = TCVMsum/nmasktot;
  PtP = Ppve'*Ppve/nmasktot; %'
  diagPtP = diag(PtP);
  pve = 100*(diagPtP.*sev_sumsqrd)/sum(diag(TCVM));
  cpve = cumsum(pve);

  fid = fopen(s.pvefile,'w');
  if(fid == -1)
    fprintf('ERROR: cannot open %s for writing\n',fname);
    return;
  end
  nn = [1:length(pve)]';%'
  fprintf(fid,'%3d   %6.4f  %6.4f\n',[nn pve cpve]'); %'
  fclose(fid);
end

%--------------------------------------------------------------%
if(~isempty(s.outvolid))
  fprintf('Starting Filtering Loop %g\n',toc);
  for sliceno = s.firstslice:s.lastslice
    fprintf('%3d ',sliceno);

    % Load the spat ev slice (again) %
    sevstem = deblank(s.sevvolid);
    [nslices nrows ncols nevs endian bext] = fmri_bvoldim(sevstem);
    if(nslices == 0)
       fprintf('ERROR: could not load %s\n',sevstem);
       return;
    end
    fname = sprintf('%s_%03d.%s',sevstem,sliceno,bext);
    sev = fmri_ldbfile(fname);
    sev = reshape(sev,[nvslice nevs]); % Dont transpose

    % Zero out the sevs to be excluded
    if(~isempty(s.sevexcl)) sev_all(:,s.sevexcl) = 0; end

    sev = sev(:,1:s.nfilter); 

    % Load the input slice %
    instem = deblank(s.involid);
    [nslices nrows ncols ntrs endian bext] = fmri_bvoldim(instem);
    if(nslices == 0)
       fprintf('ERROR: could not load %s\n',instem);
       return;
    end
    fname = sprintf('%s_%03d.%s',instem,sliceno,bext);
    yrun = fmri_ldbfile(fname);
    yrun = reshape(yrun,[nvslice ntrs])'; %'

    % Spatial Filter %
    if(~s.svestnoise) yrun2 = yrun - Pfil*sev'; %' Remove Noise Estimate
    else              yrun2 = Pfil*sev';        %' Save   Noise Estimate
    end

    yrun2b = yrun2 + repmat(yrun(1,:),[ntrs 1]); % Add mean image back
    yrun2b = reshape(yrun2b',[nrows ncols ntrs]); %'
    fname = sprintf('%s_%03d.bfloat',s.outvolid,sliceno);
    fmri_svbfile(yrun2b,fname);
    clear yrun2b;   

  end % Loop over slices
  fprintf('\n');
end

if(~isempty(s.projtcid) )
  % Pve is ntrs X nevs
  tmp = Ppve'; %'
  [n1 n2] = size(tmp);
  tmp = reshape(tmp,[1 1 n1 n2]);
  fmri_svbvolume(tmp,s.projtcid);
end

fprintf(1,'fast_evspatfil: Done %g\n',toc);

r = 0;

return;
%---\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%-----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%


%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = evspatfil_struct;
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

      case '-sev',
        arg1check(flag,narg,ninputargs);
        s.sevvolid = inputargs{narg};
        narg = narg + 1;

      case {'-nfilter','-nf'}
        arg1check(flag,narg,ninputargs);
        s.nfilter = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-sevexcl'}
        arg1check(flag,narg,ninputargs);
        s.sevexcl = [s.sevexcl sscanf(inputargs{narg},'%d',1)];
        narg = narg + 1;

      case {'-o'},
        arg1check(flag,narg,ninputargs);
        s.outvolid = inputargs{narg};
        narg = narg + 1;

      case {'-obext'},
        arg1check(flag,narg,ninputargs);
        s.obext = inputargs{narg};
        if(~strcmp(s.obext,'bshort') & ~strcmp(s.obext,'bfloat'))
  	  fprintf('ERROR: obext = %s, must be bshort or bfloat\n',...
                  s.obext);
  	  s = [];
          return;
        end
        narg = narg + 1;

      case {'-pve'},
        arg1check(flag,narg,ninputargs);
        s.pvefile = inputargs{narg};
        narg = narg + 1;

      case {'-projtc'},
        arg1check(flag,narg,ninputargs);
        s.projtcid = inputargs{narg};
        narg = narg + 1;

      case {'-polyfit'}
        arg1check(flag,narg,ninputargs);
        s.pforder = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-cutends'}
        s.cutends = 1;

      case {'-demean'}
        if(s.pforder < 0) s.pforder = 0; end
      case {'-nodemean'}
        if(s.pforder >= 0) s.pforder = -1; end
      case {'-detrend'}
        if(s.pforder < 1) s.pforder = 1; end
      case {'-nodetrend'}
        if(s.pforder >= 1) s.pforder = 0; end

      case {'-svestnoise'}
        s.svestnoise = 1; % Save noise estimate, not filtered input

      case {'-mask'}
        arg1check(flag,narg,ninputargs);
        s.maskid = inputargs{narg};
        narg = narg + 1;

      case {'-maskthresh'}
        arg1check(flag,narg,ninputargs);
        s.maskthresh = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-masksign'}
        arg1check(flag,narg,ninputargs);
        s.masksign = inputargs{narg};
        narg = narg + 1;

      case {'-automaskthresh'}
        arg1check(flag,narg,ninputargs);
        s.automaskthresh = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-automask'}
        arg1check(flag,narg,ninputargs);
        s.automaskid = inputargs{narg};
        narg = narg + 1;

      case {'-projmask'}
        arg1check(flag,narg,ninputargs);
        s.projmaskid = inputargs{narg};
        narg = narg + 1;

      case {'-projmaskthresh'}
        arg1check(flag,narg,ninputargs);
        s.projmaskthresh = sscanf(inputargs{narg},'%f',1);
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

  if(isempty(s.outvolid) & isempty(s.pvefile))
    fprintf(2,'ERROR: must spec either output vol or pve file\n');
    s=[]; return;
  end
   
  if(~isempty(s.outvolid) & isempty(s.nfilter))
    fprintf(2,'ERROR: must spec nfilter \n');
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
function s = dump_params(s,fid)
  if(nargin == 1) fid = 1; end

  fprintf(fid,'Input  Volume %s\n',s.involid);
  fprintf(fid,'SEV    Volume %s\n',s.sevvolid);

  if(~isempty(s.outvolid))
    fprintf(fid,'Output Volume %s\n',s.outvolid);
    fprintf(fid,'nfilter %d\n',s.nfilter);
  end

  if(~isempty(s.nfilter))
    fprintf(fid,'nfilter %d\n',s.nfilter);
  end

  if(~isempty(s.sevexcl))
    fprintf(fid,'excludes ');
    fprintf(fid,'%2d ',s.sevexcl);
    fprintf(fid,'\n');
  else
    fprintf(fid,'no excludes \n');
  end

  if(~isempty(s.pvefile))
    fprintf(fid,'PVE File %s\n',s.pvefile);
  end

  fprintf(fid,'PolyFit Order   %d\n',s.pforder);
  fprintf(fid,'cutends         %d\n',s.cutends);
  fprintf(fid,'firstslice   %d\n',s.firstslice);
  fprintf(fid,'nslices      %d\n',s.nslices);

return;
%--------------------------------------------------%


%------------- Print Usage ---------------------%
function print_usage(dummy)

  fprintf(1,'USAGE:\n');
  fprintf(1,'  fast_evspatfil\n');

  fprintf(1,'     -sev     sevvolid \n');
  fprintf(1,'     -i       involid  \n');

  fprintf(1,'     -o       outvolid \n');
  fprintf(1,'     -nfilter n \n');
  fprintf(1,'     -excl nth <-excl mth ...> (exclude nth sev (1-based) \n');

  fprintf(1,'     -pve pvefile \n');
  fprintf(1,'     -polyfit order  \n');

  fprintf(1,'     -mask volid \n');
  fprintf(1,'     -maskthresh thresh (0.5) \n');
  fprintf(1,'     -masksign   sign   (abs) \n');

  fprintf(1,'     -projmask volid \n');
  fprintf(1,'     -projmaskthresh thresh (0.5) \n');

return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = evspatfil_struct
  s.involid        = '';
  s.sevvolid       = '';
  s.outvolid       = '';
  s.svestnoise     = 0;
  s.obext          = '';
  s.nfilter        = [];
  s.sevexcl        = [];
  s.pvefile        = '';
  s.projtcid       = '';
  s.pforder        = -1;
  s.cutends        =  0;
  s.maskid         = '';
  s.maskthresh     = 0.5;
  s.masksign       = 'abs';
  s.mask           = [];
  s.automaskthresh = '';
  s.automaskid     = '';
  s.projmaskid     = '';
  s.projmaskthresh = 0.5;
  s.verbose = 0;
  s.firstslice = 0;
  s.nslices    = -1;
return;

