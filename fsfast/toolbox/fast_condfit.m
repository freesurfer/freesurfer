function r = fast_condfit(varargin)
% r = fast_condfit(varargin)
% Random effects fit of conditions to a weight


%
% fast_condfit.m
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

% Pre-define for compilation
hd = 0;
Nrows = 0;
Ncols = 0;
Ahat = 0;
Acov = 0;
CtC = 0;
eres_std = 0;
nsubjects = 0;

%% Parse the arguments %%
s = parse_args(varargin);
if(isempty(s)) return; end
s = check_params(s);
if(isempty(s)) return; end

ninvols = size(s.invols,1);
lastslice = s.firstslice + s.nslices - 1;

fprintf(1,'_______________________________________________\n');
fprintf(1,'Condition Fitting/Averaging Parameters\n');
condfit_print_struct(s,1);
fprintf(1,'_______________________________________________\n');

nsubjects = size(s.invols,1);
Nnnc = s.Nnnc;

if(s.synth ~= 0)
  fprintf(1,'Synthsizing with Seed %d\n',s.synth);
  randn('state',s.synth);
  rand('state', s.synth);
end

%-----------------------------------------------------------%
tic;
lastslice = s.firstslice + s.nslices - 1;

fprintf(1,'Loading Data (%d sessions)\n',nsubjects);
for slice = s.firstslice:lastslice
  fprintf('Slice %02d/%02d :',slice,s.nslices);

  havg_all = [];
  for subj = 1:nsubjects
    fprintf('%2d ',subj);
    if(rem(subj,13) == 0) fprintf('\n');end

    instem = deblank(s.invols(subj,:));

    %% Load the dat file %%
    datfile = sprintf('%s.dat',instem);
    hd = fmri_lddat3(datfile);

    %% If no condition weights are specified, use linear %%
    if(isempty(s.wcond))
      s.wcond = [1:hd.Nnnc]';  %'
      s.Nnnc = length(s.wcond);
      Nnnc = s.Nnnc;
    end

    %% Check consistency between the weights and the data file %%
    if(hd.Nnnc ~= s.Nnnc)
      msg = sprintf('Nnnc = %d, does not equal to spec = %d\n',hd.Nnnc,s.Nnnc);
      qoe(msg);error(msg);
    end

    %% Make sure excluded condition does not exceed number of conditions %%
    if(~isempty(s.exclcond))
      if(max(s.exclcond) > s.Nnnc)
        msg = sprintf('Excluded condition found out of range\n');
        qoe(msg);error(msg);
      end
    end

    % get a list of included conditions %
    tmp = ones(hd.Nnnc,1);
    tmp(s.exclcond) = 0;
    s.inclcond = find(tmp==1);
    Ninc = length(s.inclcond);

    %% Load the data %%
    fname = sprintf('%s_%03d.bfloat',instem,slice);
    havg = fast_ldsxabfile(fname);

    %hsa = fmri_ldbfile(fname);
    %% Extract averages: Nrows, Ncols,Nc,Nh  %%
    %havg  = fmri_untangle(hsa,hd);

    Nrows = size(havg,1);
    Ncols = size(havg,2);
    Nc    = size(havg,3);
    Nh    = size(havg,4);
    Nv    = Nrows * Ncols;
    Navgs = Ninc*hd.Nh;
    Nnnc  = Nc;

    %% Remove Cond0: Nrows,Ncols,Nnnc,Nh  %%
    %havg = havg(:,:,2:Nc,:);

    %%%% Synthesize Data %%%%%
    if(s.synth) 
      havg = randn(size(havg)); 
      if(1)
        havg = ones(size(havg)); 
        for nn = 1:s.Nnnc
          if(nn == 2 | nn == 4 & 0)
            havg(:,:,nn,:) = rand*havg(:,:,nn,:);
          else
            havg(:,:,nn,:) = nn*havg(:,:,nn,:);
          end
        end
      end
      hoffset = rand(Nrows,Ncols);
      %havg = havg .* repmat(hoffset,[1 1 Nnnc Nh]);
    end

    %% Percent Signal Change %%
    if(s.pctsigch)
      if(~s.synth) 
        fname = sprintf('%s-offset_%03d.bfloat',instem,slice);
        hoffset = fmri_ldbfile(fname);
      end
      ind = find(hoffset==0);
      hoffset(ind) = 10^10;
      havg = 100*havg ./ repmat(hoffset,[1 1 s.Nnnc hd.Nh]);
    end

    %% Zero Prestim %%
    if(s.zprestim)
      nprestim = round(hd.TPreStim/hd.TR) + 1;
      hprestim = mean(havg(:,:,:,1:nprestim),4);
      havg = havg - repmat(hprestim,[1 1 1 hd.Nh]);
    end

    %% Collect vertices to save %%
    if(~isempty(s.vtxlist))
      if(subj==1) 
        vtxsave = zeros(s.Nvl,nsubjects,Nnnc,Nh);
        [vtxsavecol vtxsaverow] = ind2sub([Ncols Nrows],s.vtxlist+1);
       end
      tmp = permute(havg, [2 1 3 4]);
      tmp = reshape(tmp, [Nv Nnnc Nh]);
      vtxsave(:,subj,:,:) = tmp(s.vtxlist+1,:,:);
      clear tmp;
      %vtxsave(:,subj,:,:) = havg(vtxsaverow,vtxsavecol,:,:);
    end

    % Only keep the ones that should be included %
    havg = havg(:,:,s.inclcond,:);

    % Keep a list of all the data %
    % havg_all = [havg_all; havg];
    havg_all(subj,:,:,:,:) = havg;
  end
  fprintf('  %f\n',toc);

  %% Reshape data for deconvolution %
  % Starts off: Nj Nr   Nc Ninc Nh
  % Permute to: Nj Ninc Nr Nc   Nh 
  havg_all = permute(havg_all, [1 4 2 3 5]);
  % Reshape to: Nj*Ninc Nr*Nc*Nh 
  havg_all = reshape(havg_all, [nsubjects*Ninc Nrows*Ncols*hd.Nh]);

  %% Create the convolution matrix %%
  cc = reshape1d(repmat(s.inclcond',[nsubjects 1]));  %'
  w = s.wcond(cc);
  if(s.fitorder == 1)  C = [ones(size(w)) w ];
  else                 C = [ones(size(w)) w w.^2];
  end
  Nfit = s.fitorder+1;

  %% Create the de-convolution matrix %%
  CtC = C'*C; %'
  fprintf(1,'Condition of deconv matrix: %g\n',cond(CtC));

  Acov = inv(CtC);

  %% Deconvolve %%
  Ahat = Acov*C'*havg_all; %'

  %% Compute the Residual Error and its variance %%
  Hhat = C*Ahat;
  eres = Hhat - havg_all;
  eres_std = std(eres);

  %% Dont divide by 0 %%
  ind = find(eres_std == 0);
  eres_std(ind) = 10^10;

  %% Compute t-values %%
  dof = nsubjects*Ninc-s.fitorder-1;
  d = repmat(diag(Acov),[1 Nv*hd.Nh']); %'
  t = (Ahat./repmat(eres_std,[size(Ahat,1) 1]))./sqrt(d);
  sig = tTest(dof,reshape1d(t));
  sig = reshape(sig, size(t));

  %% Reshape data back to human-readable form %%
  % Reshape to: Nfit Nr Nc Nh 
  sig = reshape(sig, [Nfit Nrows Ncols Nh]);
  % Permute to: Nr Nc Nh Nfit
  sig = permute(sig, [2 3 4 1]);
  % Get rid of offset
  sig = sig(:,:,:,2);

  %% Save as log10 %%
  sig = sign(sig) .* abs(log10(sig));
  outfile = sprintf('%s_%03d.bfloat',s.outvol,slice);
  fmri_svbfile(sig, outfile);

  if(~isempty(s.eresstdvol))
    fprintf('Saving residual standard deviation\n');
    % Reshape %
    eres_std  = reshape(eres_std, [1 Nrows Ncols Nh]);
    eres_std  = permute(eres_std, [2 3 4 1]);
    eres_std0 = eres_std; % keep the original one around
    outfile = sprintf('%s_%03d.bfloat',s.eresstdvol,slice);
    fmri_svbfile(eres_std, outfile);
  end

  if(~isempty(s.coeffvol))
    fprintf('Saving linear coefficients\n');
    % Reshape %
    Ahat  = reshape(Ahat, [Nfit Nrows Ncols Nh]);
    Ahat  = permute(Ahat, [2 3 4 1]);
    Ahat0 = Ahat; % keep the original one around
    Ahat  = Ahat(:,:,:,2);
    outfile = sprintf('%s_%03d.bfloat',s.coeffvol,slice);
    fmri_svbfile(Ahat, outfile);
  end

  % Save vertex list %
  if(~isempty(s.vtxlist))
    fprintf(1,'Saving %d vertices\n',s.Nvl);

    % save summary file %
    fname = sprintf('%s.sum',s.vtxfile);
    fidsum = fopen(fname,'w');    

    % get coefficients %
    tmp = permute(Ahat0,[2 1 3 4]);
    tmp = reshape(tmp, [Nv Nh Nfit]);
    vtxcoef = tmp(s.vtxlist+1,:,:);

    % get signficances %
    tmp = permute(sig,[2 1 3]);
    tmp = reshape(tmp, [Nv Nh]);
    vtxsig = tmp(s.vtxlist+1,:);

    for vtx = 1:s.Nvl
      fname = sprintf('%s-%d.dat',s.vtxfile,s.vtxlist(vtx));
      fid = fopen(fname,'w');
      for subj = 1:nsubjects,
        for c = 1:Nnnc,
          fprintf(fid,'%2d %2d %4.2f ',subj,c,s.wcond(c));
          fprintf(fid,'%7.4f ',vtxsave(vtx,subj,c,:));
          fprintf(fid,'\n');
        end
      end
      fclose(fid);
 
      fprintf(fidsum,'%5d ',s.vtxlist(vtx));
      for h = 1:Nh
        fprintf(fidsum,'%7.4f %7.4f %7.4f ',...
	      vtxcoef(vtx,h,1),vtxcoef(vtx,h,2),vtxsig(vtx,h));
      end % loop over Nh
      fprintf(fidsum,'\n');

    end % loop over vtx

  end

end % loop over slices 

fname = sprintf('%s.cflog',s.outvol);
fid = fopen(fname,'w');
condfit_print_struct(s,fid);
fclose(fid);

fprintf('Elapsed Time: %g\n',toc);

r = 0;
fprintf(1,'fast_condfit: completed SUCCESSFULLY\n');

return;
%----------------------------------------------------------%
%----------------------------------------------------------%


%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = condfit_struct;
  inputargs = varargin{1};
  ninputargs = length(inputargs);

  narg = 1;
  while(narg <= ninputargs)

    flag = deblank(inputargs{narg});
    narg = narg + 1;
    %fprintf(1,'Argument: %2d %s\n',narg,flag);
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

      case {'-sig','-o'},
        arg1check(flag,narg,ninputargs);
        s.outvol = inputargs{narg};
        narg = narg + 1;

      case {'-eresstd'}
        arg1check(flag,narg,ninputargs);
        s.eresstdvol = inputargs{narg};
        narg = narg + 1;

      case {'-coeff'}
        arg1check(flag,narg,ninputargs);
        s.coeffvol = inputargs{narg};
        narg = narg + 1;

      case '-fitorder',
        arg1check(flag,narg,ninputargs);
        s.fitorder = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-wcond',
        arg1check(flag,narg,ninputargs);
        s.Nnnc = 0;
        while(narg <= ninputargs & ~strncmp(inputargs{narg},'-',1) )
          s.Nnnc = s.Nnnc + 1;
          s.wcond(s.Nnnc) = sscanf(inputargs{narg},'%f',1);
          narg = narg + 1;
        end
        if(s.Nnnc == 0) 
          fprintf(1,'ERROR: no condition weights listed\n');
          s = [];
          return;
        end
        s.wcond = s.wcond';%'

      case '-vtxlist',
        arg1check(flag,narg,ninputargs);
        s.Nvl = 0;
        while(narg <= ninputargs & ~strncmp(inputargs{narg},'-',1) )
          s.Nvl = s.Nvl + 1;
          s.vtxlist(s.Nvl) = sscanf(inputargs{narg},'%f',1);
          narg = narg + 1;
        end
        if(s.Nvl == 0) 
          fprintf(1,'ERROR: no vertices listed\n');
          s = [];
          return;
        end

      case '-vtxfile',
        arg1check(flag,narg,ninputargs);
        s.vtxfile = inputargs{narg};
        narg = narg + 1;

      case '-vtxframe',
        arg1check(flag,narg,ninputargs);
        s.vtxframe = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-exclcond',
        arg1check(flag,narg,ninputargs);
        s.Nexcl = 0;
        while(narg <= ninputargs & ~strncmp(inputargs{narg},'-',1)  )
          s.Nexcl = s.Nexcl + 1;
          s.exclcond(s.Nexcl) = sscanf(inputargs{narg},'%f',1);
          narg = narg + 1;
        end
        if(s.Nexcl == 0) 
          fprintf(1,'ERROR: no excluded conditions listed\n');
          s = [];
          return;
        end

      case '-logfile',
        arg1check(flag,narg,ninputargs);
        s.logfile = inputargs{narg};
        narg = narg + 1;

      case '-pctsigch', % 
        s.pctsigch = 1;

      case '-zprestim', % 
        s.zprestim = 1;

      case '-monly', % ignore
        arg1check(flag,narg,ninputargs);
        s.mfile = inputargs{narg};
        narg = narg + 1;

      case '-synth',
        arg1check(flag,narg,ninputargs);
        s.synth = sscanf(inputargs{narg},'%d',1);
        if(s.synth == -1) s.synth = sum(100*clock); end
        narg = narg + 1;

      case '-nolog',
        s.nolog = 1;
        s.logfid = 0;

      case '-debug', % ignore

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

  if(isempty(s.outvol))
    fprintf(1,'ERROR: no output specified\n');
    error;
  end

  if(~isempty(s.vtxlist) & isempty(s.vtxfile))
    fprintf(1,'ERROR: vtxfile must be specified with vtxlist\n');
    error;    
  end

return;

%--------------------------------------------------%
%% Print Usage 
function print_usage
  fprintf(1,'USAGE:\n');
  fprintf(1,'  fast_condfit\n');
  fprintf(1,'     -i invol1 -i invol2 ... \n');
  fprintf(1,'     -sig        volid\n');
  fprintf(1,'     -eresstd    volid \n');
  fprintf(1,'     -coeff      volid \n');
  fprintf(1,'     -firstslice slice\n');
  fprintf(1,'     -nslices    nslices\n');
  fprintf(1,'     -fitorder   order (<1>,2)\n');
  fprintf(1,'     -wcond      weightlist \n');
  fprintf(1,'     -exclcond   condition list \n');
  fprintf(1,'     -pctsigch   use percent signal change  \n');
  fprintf(1,'     -zprestim   zero prestimulus window  \n');
  fprintf(1,'     -synth  \n');
return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = condfit_struct
  s.invols     = '';
  s.involfmt   = '';
  s.involtype  = '';
  s.firstslice = 0;
  s.nslices    = -1;
  s.outvol     = '';
  s.eresstdvol = '';
  s.coeffvol   = '';
  s.wcond      = [];
  s.Nnnc       = 0;
  s.Nexcl      = 0;
  s.exclcond   = [];
  s.Nframes    = 0;
  s.fitorder   = 1;
  s.logfile    = '';
  s.logfid     = 0;
  s.nolog      = 0;
  s.synth      = 0;
  s.pctsigch   = 0; % fit percent signal change
  s.zprestim   = 0; % zero prestim baseline 
  s.Nvl        = 0;  % number of vertices in list
  s.vtxlist    = []; % list of vertices to save in text file
  s.vtxfile    = ''; % file name
  s.vtxframe   = 0;  % frame to save 
return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Print data structure
function s = condfit_print_struct(s,fid)
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
  fprintf(fid,'outvol      %s\n',s.outvol);

  if(~isempty(s.eresstdvol))
    fprintf(fid,'eresstd      %s\n',s.eresstdvol);
  end

  if(~isempty(s.coeffvol))
    fprintf(fid,'coeff      %s\n',s.coeffvol);
  end

  fprintf(fid,'wcond         ');
  fprintf(fid,'%f ',s.wcond);
  fprintf(fid,'\n');

  fprintf(fid,'excluded conditions ');
  fprintf(fid,'%d ',s.exclcond);
  fprintf(fid,'\n');

  fprintf(fid,'pctsigch     %d\n',s.pctsigch);
  fprintf(fid,'zprestim     %d\n',s.zprestim);

  fprintf(fid,'fitorder     %d\n',s.fitorder);
  fprintf(fid,'logfile     %s\n',s.logfile);
  fprintf(fid,'nolog       %d\n',s.nolog);

  fprintf(fid,'vtxlist ');
  fprintf(fid,'%d ',s.vtxlist);
  fprintf(fid,'\n');
  fprintf(fid,'vtxfile   %s\n',s.vtxfile);
  fprintf(fid,'vtxframe  %d\n',s.vtxframe);

  fprintf(fid,'synth       %d\n',s.synth);

return;
%--------------------------------------------------%




