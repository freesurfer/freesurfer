function r = fast_func2roi(varargin)
% r = fast_func2roi(varargin)


%
% fast_func2roi.m
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

version = 'fast_func2roi.m @FS_VERSION@';
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

if(s.verbose)
  f2r_print_struct(s);
end

%%%%  Get info about the input functional volume %%%%
[nslices nrows ncols nt endian bext hdrdat] = fmri_bvoldim(s.funcvolid);

% Load the label2func registration matrix %
[s.l2f subject inres betres ireg] = fmri_readreg(s.l2f_regfile);

%------- Contruct the functional quantization matrix -------%
fvox = [inres inres betres];
fdim = [nrows ncols nslices];
Qf = zeros(4,4);
Qf = [-fvox(1)      0       0     fvox(1)*(fdim(1)-1)/2; ...
         0          0    fvox(3) -fvox(3)*(fdim(3)-1)/2; ...
         0      -fvox(2)    0     fvox(2)*(fdim(2)-1)/2; ...
         0          0       0       1];

%% ------------  Load the label data -------------------- %
[vtxno xyz vstat msg] = fast_ldlabel(s.labelfile);
if( ~strcmp(msg,'OK') )
  qoe(msg); error(msg);
end

% --- Compute the functional subscripts of the label, 0-based -----%
rcs_f = (inv(Qf)*s.l2f*([xyz ones(size(xyz,1),1)])')';
oob = find( rcs_f(:,1) < 0 | rcs_f(:,1) > (nrows-1) | ...
            rcs_f(:,2) < 0 | rcs_f(:,2) > (ncols-1) | ...
            rcs_f(:,3) < 0 | rcs_f(:,3) > (nslices-1) );
ntmp = ones(size(xyz,1),1);
ntmp(oob) = 1;
itmp = find(ntmp==1);
rcs_f= rcs_f(itmp,1:3);

% --- Get the indices and weights of the label for trilinear interp --- %
% Note conversion to one-based %
[ind wind] = tliweights([nrows ncols nslices], rcs_f + 1);
ind  = reshape1d(ind);
wind = reshape1d(wind);

%--------------------------------------------------------- %
% --------------------  Mask ----------------------------- % 
if(~isempty(s.maskvolid))

  % -------------  Load the mask ---------------------- % 
  fprintf(1,'Loading mask volume\n');
  [maskfv msg] = fast_ldmri(s.maskvolid,[],[],[],[],s.mplane+1);
  if( ~strcmp(msg,'OK') ) qoe(msg); error(msg); end

  % Get indices of suprathreshold voxels in mask volume %
  switch(s.msign)
     case {'abs'}
       ind_mm = find(abs(maskfv.data) > s.mthresh);
     case {'pos'}
       ind_mm = find(maskfv.data > s.mthresh);
     case {'neg'}
       ind_mm = find(-maskfv.data > s.mthresh);
  end

  % Make sure that there are voxels in the mask %
  if(isempty(ind_mm))
    msg = sprintf('No (%s) voxels found above threshold %f',s.msign,s.mthresh);
    qoe(msg);error(msg);
  end
  NMaskVoxels = length(ind_mm);

  % Get the subscripts of the indices (one-based) %
  [r_mm c_mm s_mm] = ind2sub(size(maskfv.data),ind_mm);

  % Convert to zero-based %
  r_mm = r_mm - 1;
  c_mm = c_mm - 1;
  s_mm = s_mm - 1;

  % Assign the mask-to-functional registration %
  if(~isempty(s.m2f_regfile))
    % Load the registration file %
    [s.m2f m_subject m_inres m_betres m_ireg] = fmri_readreg(s.m2f_regfile);

    % Make sure it is consistent with the other register.dat %
    if(m_inres ~= inres | m_betres ~= betres)
      msg = sprintf('Mask RegFile inconsistent with Label RegFile');
      qoe(msg); error(msg);
    end

  else
    % No registration file, assume identity %
    s.m2f = eye(4);    
  end

  % If left undefined, set the size of the mask voxel to that
  % of the functional voxel.
  if(isempty(s.mvoxsize)) s.mvoxsize = [inres inres betres]; end

  % Create the mask volume quanization matrix %
  mvox = s.mvoxsize;
  mdim = size(maskfv.data);
  Qm = [-mvox(1)      0       0     mvox(1)*(mdim(1)-1)/2; ...
           0          0    mvox(3) -mvox(3)*(mdim(3)-1)/2; ...
           0      -mvox(2)    0     mvox(2)*(mdim(2)-1)/2; ...
           0          0       0       1];

  % Compute the location of the mask voxels in functional subscript space %
  % Zero-based
  rcs_mf = (inv(Qf)*s.m2f*Qm*([r_mm c_mm s_mm ones(NMaskVoxels,1)])')';
  rcs_mf = round(rcs_mf(:,1:3));

  % Remove the mask voxels that are outsize the functional volume
  oob = find( rcs_mf(:,1) < 0 | rcs_mf(:,1) > (nrows-1) | ...
              rcs_mf(:,2) < 0 | rcs_mf(:,2) > (ncols-1) | ...
              rcs_mf(:,3) < 0 | rcs_mf(:,3) > (nslices-1) );
  ntmp = ones(size(xyz,1),1);
  ntmp(oob) = 1;
  itmp = find(ntmp==1);

  % Make sure that there are voxels left in the mask %
  if(isempty(itmp))
    msg = sprintf('No (%s) voxels found above threshold %f',s.msign,s.mthresh);
    qoe(msg);error(msg);
  end
  rcs_mf= rcs_mf(itmp,:);

  % Compute the functional indices of the mask. %
  % Note conversion to one-based.
  ind_mf = sub2ind([nrows ncols nslices], rcs_mf(:,1)+1, ...
                   rcs_mf(:,2)+1, rcs_mf(:,3)+1);

  % Just keep the uniqe ones %
  ind_mf = unique(ind_mf);

  fprintf(1,'INFO: found %d mask voxels\n',length(ind_mf));

  % Get intersection between label and mask %
  [tmp ind_ind tmp2] = intersect(ind,ind_mf);
  tmpind = ind;
  ind  = ind(ind_ind);
  wind = wind(ind_ind);

end

nunique = length(unique(ind));
fprintf(1,'INFO: There are %d unique functional voxels in ROI\n',nunique);

fprintf(1,'INFO: Loading functional volume\n');
[funcfv msg] = fast_ldmri(s.funcvolid);
if( ~strcmp(msg,'OK') ) qoe(msg); error(msg); end
nv = prod(funcfv.szvol(1:3));

ftmp = reshape(funcfv.data, [nv funcfv.szvol(4)]);
yroi = mean(ftmp(ind,:));


r = 0;
fprintf(1,'fast_func2roi: completed SUCCESSFULLY\n');

keyboard

return;
%----------------------------------------------------------%
%----------------------------------------------------------%
%----------------------------------------------------------%
%----------------------------------------------------------%


%--------------------------------------------------%
%% Print Usage 
function print_usage(dummy)

  fprintf(1,'USAGE:\n');
  fprintf(1,'  fast_func2roi\n');
  fprintf(1,'     -f        funcvolid \n');
  fprintf(1,'     -l        labelfile \n');
  fprintf(1,'     -l2f      reg-label2func \n');
  fprintf(1,'     -m        maskvolid \n');
  fprintf(1,'     -m2f      reg-mask2func  \n');
  fprintf(1,'     -mvoxsize row_sz col_sz slice_sz \n');
  fprintf(1,'     -mthresh  masktreshold \n');
  fprintf(1,'     -msign    masksign (abs,pos,neg)\n');
  fprintf(1,'     -mplane   maskplane (zero-based)\n');
  fprintf(1,'     -o        outvolid \n');

return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = f2r_struct
  s.funcvolid   = '';

  s.labelfile   = '';
  s.l2f_regfile = '';
  s.l2f         = [];

  s.maskvolid   = '';
  s.m2f_regfile = '';
  s.m2f         = [];
  s.mvoxsize    = [];
  s.mthresh     = 0.5;
  s.msign       = 'abs';
  s.mplane      = 0;

  s.outvolid   = '';

  s.verbose = 0;
return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Print data structure
function s = f2r_print_struct(s,fid)
  if(nargin == 1) fid = 1; end

  fprintf(fid,'funcvolid    %s\n',s.funcvolid);
  fprintf(fid,'labelfile    %s\n',s.labelfile);
  fprintf(fid,'l2f_regfile  %s\n',s.l2f_regfile);

  fprintf(fid,'maskvolid    %s\n',s.maskvolid);
  fprintf(fid,'m2f_regfile  %s\n',s.m2f_regfile);
  fprintf(fid,'maskvoxsize  %f %f %f\n',s.mvoxsize);
  fprintf(fid,'mask thresh  %f\n',s.mthresh);
  fprintf(fid,'mask sign    %s\n',s.msign);
  fprintf(fid,'mask plane   %d\n',s.mplane);

  fprintf(fid,'outvolid     %s\n',s.outvolid);

return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Check argument consistency, etc %%%
function s = check_params(s)

  fprintf(1,'Checking Parameters\n');

  if(size(s.funcvolid,1) < 1) 
    msg = sprintf('ERROR: must have an input functional volume');
    qoe(msg);error(msg);
  end

  if(size(s.labelfile,1) < 1) 
    msg = sprintf('ERROR: must specify a label file');
    qoe(msg);error(msg);
  end

  if(size(s.l2f_regfile,1) < 1) 
    msg = sprintf('ERROR: must specify a label register file');
    qoe(msg);error(msg);
  end

  if(~ (strcmp(s.msign,'abs') | strcmp(s.msign,'pos') | strcmp(s.msign,'neg')))
    msg = sprintf('ERROR: -msign must be either abs, pos, or neg');
    qoe(msg);error(msg);
  end


return;



%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = f2r_struct;
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

      case '-f',
        arg1check(flag,narg,ninputargs);
        s.funcvolid = inputargs{narg};
        narg = narg + 1;

      case {'-l'}
        arg1check(flag,narg,ninputargs);
        s.labelfile = inputargs{narg};
        narg = narg + 1;

      case {'-l2f'}
        arg1check(flag,narg,ninputargs);
        s.l2f_regfile = inputargs{narg};
        narg = narg + 1;

      case '-m',
        arg1check(flag,narg,ninputargs);
        s.maskvolid = inputargs{narg};
        narg = narg + 1;

      case {'-m2f'}
        arg1check(flag,narg,ninputargs);
        s.m2f_regfile = inputargs{narg};
        narg = narg + 1;

      case {'-mthresh'}
        arg1check(flag,narg,ninputargs);
        s.mthresh = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-mvoxsize'}
        arg3check(flag,narg,ninputargs);
        s.mvoxsize = sscanf(inputargs{narg},'%f',3);
        narg = narg + 3;

      case {'-msign'}
        arg1check(flag,narg,ninputargs);
        s.msign = inputargs{narg}; 
        narg = narg + 1;

      case {'-mplane'}
        arg1check(flag,narg,ninputargs);
        s.mplane = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-o',
        arg1check(flag,narg,ninputargs);
        s.outvolid = inputargs{narg};
        narg = narg + 1;

      case '-monly', % ignore
        arg1check(flag,narg,ninputargs);
        narg = narg + 1;

      case '-verbose',
        s.verbose = 1;

      case {'-debug','-echo','-umask'}, % ignore

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
%% Check that there are at least three more arguments %%
function arg3check(flag,nflag,nmax)
  if(nflag+2 > nmax) 
    fprintf(1,'ERROR: Flag %s needs one argument',flag);
    error;
  end
return;

