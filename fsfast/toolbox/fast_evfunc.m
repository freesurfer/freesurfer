function r = fast_evfunc(varargin)
% r = fast_evfunc(varargin)
%
% Computes the spatial eigen structure for a list of input 
% functional volumes. If a mask is supplied (or automask is 
% specified, the input volumes will be masked off. After
% masking, temporal covmtx is computed and saved. The Temporal 
% EVects and EVals are then computed from the tcvm. These
% are then used to compute the Spatial EVects. The PVE and CPVE 
% of the  EVs are also computed and saved in outstem-pve.dat.
%
% Masking. If a mask is specified, all the voxels over 0.5
% or under -0.5 are used. The threshold can be changed as can
% the sign.
%
% Automasking. If -automaskthresh T is specified, only voxels whose
% temporal mean is greater than T times the global mean will be 
% used. The resulting mask can also be saved.
%
%


%
% fast_evfunc.m
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

version = 'fast_evfunc.m @FS_VERSION@';
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

main_print_struct(s);

% Get basic info from the first run %
instem = deblank(s.invollist(1,:));
[nslices nrows ncols ntrs] = fmri_bvoldim(instem);
nv = nslices*nrows*ncols;
nvslice = nrows*ncols;

if(nv == 0)
  fprintf('ERROR: reading %s\n',instem);
  return;
end

tic;

%----------------- masking ------------------------------%
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
  instem = deblank(s.invollist(1,:));
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
    fmri_svbvolume(s.mask,s.automaskid,[],'bshort');
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
%----------------- done masking ----------------------------------%

%--------------------------------------------------------------%
% Create Temporal Smoothing Matrix %
fprintf('tsmooth = %d\n',s.tsmooth);
if(s.tsmooth > 0)
  fprintf('Creating Temporal Smoothing Matrix\n');
  TSmoothMtx = zeros(ntrs/2,ntrs);
  TSmoothMtx(:,[1 2]) = 1;
  for n = 2:2:(ntrs-2)
    rw = n/2 + 1;
    TSmoothMtx(rw,:) = fast_mshift(TSmoothMtx(rw,:),[0 n],0);
  end
  ntrs = ntrs/2;
else
  TSmoothMtx = [];
end

% ------ Construct DeMeaning/DeTrending Matrix ----------- %
if(s.pforder > 0)
  X = fast_polytrendmtx(1,ntrs,1,s.pforder);
  D = eye(ntrs) - X*inv(X'*X)*X';
  nkeep = ntrs - (s.pforder+1);
else
  D = [];
  nkeep = ntrs;
end

if( ~isempty(s.nkeep) )
  if(nkeep > s.nkeep) nkeep = s.nkeep; end 
end
fprintf('nkeep = %d\n',nkeep);

%--------------------------------------------------------------%
fprintf('Computing Temporal Covariance Matrix %g\n',toc);
TCVMsum = 0;
for sliceno = s.firstslice:s.lastslice
  fprintf('%3d ',sliceno);

    % Mask indicies for a given slices %
    imask = find(s.mask(sliceno+1,:,:)==1);
    nmask = length(imask);

    if(nmask > 0)
      yall = [];
      for runno = 1:s.nruns
        instem = deblank(s.invollist(runno,:));
        [nslices nrows ncols ntrs endian bext] = fmri_bvoldim(instem);
        if(nslices == 0)
          fprintf('ERROR: could not load %s\n',instem);
          return;
        end
        fname = sprintf('%s_%03d.%s',instem,sliceno,bext);
        yrun = fmri_ldbfile(fname);
        yrun = reshape(yrun,[nvslice ntrs])'; %'
        yrun = yrun(:,imask);
        if(~isempty(TSmoothMtx)) yrun = TSmoothMtx*yrun; end
        if(~isempty(D)) yrun = D*yrun; end % Remove drift
        yall = [yall; yrun];
      end % Loop over runs %
    
      TCVMslice = (yall*yall'); %' Temp Cov Mtx
      TCVMsum = TCVMsum + TCVMslice; % Accumulate over all slices

    end %if(nmask > 0)%

end % Loop over slices
fprintf('\n');
TCVM = TCVMsum/nmasktot;

% Save TCVM %
fname = sprintf('%s-tcvm.bfloat',s.outstem);
fmri_svbfile(TCVM,fname);
fprintf('TCVM Done %g\n',toc);

fprintf('Computing Temporal SVD\n');
[TEigVects EigVals tmp] = svd(TCVM);

% Compute Percent Variance Explained by each EV
pve = 100*diag(EigVals)/sum(diag(EigVals));
cpve = cumsum(pve); % cumulative
fname = sprintf('%s-pve.dat',s.outstem);
fid = fopen(fname,'w');
if(fid == -1)
  fprintf('ERROR: cannot open %s for writing\n',fname);
  return;
end
nn = [1:length(pve)]';%'
fprintf(fid,'%3d   %6.4f   %6.4f\n',[nn pve cpve]'); %'
fclose(fid);

% Do not keep all of them (some will be zero) %
TEigVects = TEigVects(:,1:nkeep);
EigVals  = EigVals(1:nkeep,1:nkeep);

%--------------------------------------------------------------%
fprintf('Computing Spatial Eigen Vectors %g\n',toc);
SEigVectsSum2 = 0;
for sliceno = s.firstslice:s.lastslice
  fprintf('%3d ',sliceno);

  % Mask indicies for a given slices %
  imask = find(s.mask(sliceno+1,:,:)==1);
  nmask = length(imask);

  if(nmask > 0)
    yall = [];
    for runno = 1:s.nruns
      instem = deblank(s.invollist(runno,:));
      [nslices nrows ncols ntrs endian bext] = fmri_bvoldim(instem);
      if(nslices == 0)
        fprintf('ERROR: could not load %s\n',instem);
        return;
      end
      fname = sprintf('%s_%03d.%s',instem,sliceno,bext);
      yrun = fmri_ldbfile(fname);
      yrun = reshape(yrun,[nvslice ntrs])'; %'
      yrun = yrun(:,imask);
      if(~isempty(TSmoothMtx)) yrun = TSmoothMtx*yrun; end
      if(~isempty(D)) yrun = D*yrun; end % Remove mean/trend 
      yall = [yall; yrun];
    end % Loop over runs %
    
    % Compute spatial eigenvectors %
    SEigVects = yall'*(TEigVects*sqrt(inv(EigVals*nmasktot))); %'
    SEigVectsSum2 = SEigVectsSum2 + sum(SEigVects.^2);

    % Insert zeros for voxels not in the mask %
    tmp = zeros(nvslice,nkeep);
    tmp(imask,:) = SEigVects;
    SEigVects = tmp;
    clear tmp;
    
    % Reshape so that it can be saved %
    SEigVects = reshape(SEigVects, [nrows ncols nkeep]);

    % Save %
    fname = sprintf('%s_%03d.bfloat',s.outstem,sliceno);
    fmri_svbfile(SEigVects,fname);

  else % Create a blank slice as filler
    SEigVects = zeros(nrows,ncols,nkeep);
    fname = sprintf('%s_%03d.bfloat',s.outstem,sliceno);
    fmri_svbfile(SEigVects,fname);
  end %if(nmask > 0)%

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
  fprintf(1,'  fast_evfunc\n');
  fprintf(1,'     -i invol ... \n');
  fprintf(1,'     -o outstem : output will be outstem, -pve, -tcvm \n');
  fprintf(1,'     -nkeep n : number of eigenvectors to keep (default all) \n');
  fprintf(1,'     -polyfit order : fit nth order polynomial \n');

  fprintf(1,'     -mask volid \n');
  fprintf(1,'     -maskthresh thresh (0.5) \n');
  fprintf(1,'     -masksign   sign   (abs) \n');

  fprintf(1,'     -automaskthresh rthresh \n');
  fprintf(1,'     -automask volid : save automask to volid\n');
return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = evfunc_struct
  s.invollist      = '';
  s.nruns          = 0;
  s.tsmooth        = 0;
  s.outstem        = '';
  s.nkeep          = [];
  s.pforder        = -1;
  s.maskid         = '';
  s.maskthresh     = 0.5;
  s.masksign       = 'abs';
  s.automaskthresh = '';
  s.automaskid     = '';
  s.mask           = [];
  s.verbose = 0;
  s.cutends = 0;
  s.firstslice = 0;
  s.nslices    = -1;
return;

%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = evfunc_struct;
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

      case {'-o'},
        arg1check(flag,narg,ninputargs);
        s.outstem = inputargs{narg};
        narg = narg + 1;

      case {'-nk','-nkeep'}
        arg1check(flag,narg,ninputargs);
        s.nkeep = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-polyfit'}
        arg1check(flag,narg,ninputargs);
        s.pforder = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-demean'}
        if(s.pforder < 0) s.pforder = 0; end
      case {'-nodemean'}
        if(s.pforder >= 0) s.pforder = -1; end
      case {'-detrend'}
        if(s.pforder < 1) s.pforder = 1; end
      case {'-nodetrend'}
        if(s.pforder >= 1) s.pforder = 0; end

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

      case '-cutends',
        s.cutends = 1;

      case {'-fs','-firstslice'}
        arg1check(flag,narg,ninputargs);
        s.firstslice = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-ns','-nslices'}
        arg1check(flag,narg,ninputargs);
        s.nslices = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-tsmooth'}
        arg1check(flag,narg,ninputargs);
        s.tsmooth = sscanf(inputargs{narg},'%d',1);
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

  if(s.nruns < 1) 
    fprintf(2,'ERROR: No input volumes specified\n');
    s=[]; return;
  end

  if(isempty(s.outstem))
    fprintf(2,'ERROR: No outstem specified\n');
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
function s = main_print_struct(s,fid)
  if(nargin == 1) fid = 1; end

  fprintf(fid,'Number of Runs: %d\n',s.nruns);
  fprintf(fid,'Input Volume List\n');
  for n = 1:size(s.invollist,1),
    fprintf(fid,'  %d  %s\n',n,s.invollist(n,:));    
  end

  fprintf(fid,'Outstem %s\n',s.outstem);
  fprintf(fid,'poly order   %d\n',s.pforder);
  fprintf(fid,'cutends      %d\n',s.cutends);
  fprintf(fid,'firstslice   %d\n',s.firstslice);
  fprintf(fid,'nslices      %d\n',s.nslices);
  fprintf(fid,'tsmooth      %d\n',s.tsmooth);

return;
%--------------------------------------------------%


