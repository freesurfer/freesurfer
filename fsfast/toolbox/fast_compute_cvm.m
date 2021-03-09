function r = fast_compute_cvm(varargin)
% r = fast_compute_cvm(varargin)
%
% Computes the temporal covariance matrix (CVM) of a functional volume.
% Options:
%  - compute CVM from a set of voxels specified in mask
%  - compute complementary CVM from the set of voxels not 
%    specified in mask
%  - specify the threhsold for mask
%  - remove mean  before computing CVM (default)
%  - remove trend before computing CVM (default)
%
% Saves results in CVM format (.bfloat, .hdr, and .cvm files)
%
% See also: fmri_cvmstruct, fmri_svcvm, fmri_ldcvm
%  


%
% fast_compute_cvm.m
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

version = 'fast_compute_cvm.m @FS_VERSION@';
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

fprintf('Loading data\n');
f = fmri_ldbvolume(s.involid);
if(isempty(f)) return; end
[nslices nrows ncols ntrs] = size(f);
nv = nslices*nrows*ncols;
f = reshape(f, [nv ntrs])';%'

if(s.synth)
  fprintf('Synthsizing data\n');
  f = randn(size(f));
end

% Load mask %
if(~isempty(s.maskid))  
  fprintf('Loading mask\n');
  mask = fmri_ldbvolume(s.maskid); 
  if(isempty(mask)) return; end
  iover  = find(mask >  s.thresh);
  iunder = find(mask <= s.thresh);
else
  mask = [];
  iover = [1:nv];
  iunder = [];
end

nover  = length(iover);
nunder  = length(iunder);
fprintf('nover = %6d, nunder = %6d\n',nover,nunder);

if(~isempty(s.ccvmstem) & nunder == 0)
  fprintf('ERROR: cannot compute complementary cvm.\n');
  fprintf('       No voxels under threshold.\n');
  return;
end

% Remove mean and/or trend %
if(~ s.keepmean | ~ s.keeptrend)
  fprintf('Removing mean (%d) and/or trend (%d)\n',...
          ~s.keepmean,~s.keeptrend)
  X = [];
  Xbaseline = fast_baselinemtx(1,ntrs,1);
  Xtrend    = fast_trendmtx(1,ntrs,1);
  if(~ s.keepmean)  X = [X Xbaseline]; end
  if(~ s.keeptrend) X = [X Xtrend];    end
  bhat = inv(X'*X)*X'*f;
  f = f - X*bhat;
end

fprintf('Computing CVM \n');
cvmover = fmri_cvmstruct;
cvmover.cvm = f(:,iover) * f(:,iover)'/nover; %'
cvmover.n = nover;
fmri_svcvm(cvmover,s.cvmstem);

if(~isempty(s.ccvmstem))
  fprintf('Computing Complementary CVM \n');
  cvmunder = fmri_cvmstruct;
  cvmunder.cvm = f(:,iunder) * f(:,iunder)'/nunder; %'
  cvmunder.n = nunder;
  fmri_svcvm(cvmunder,s.ccvmstem);
end

return;
%---\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%-----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%


%------------- Print Usage ---------------------%
function print_usage(dummy)

  fprintf(1,'USAGE:\n');
  fprintf(1,'  fast_compute_cvm\n');
  fprintf(1,'     -i   invol  \n');
  fprintf(1,'     -cvm cvmstem \n');
  fprintf(1,'     -ccvm complementary cvmstem  \n');
  fprintf(1,'     -mask  maskvol \n');
  fprintf(1,'     -thresh mask threshold \n');
  fprintf(1,'     -keepmean \n');
  fprintf(1,'     -keeptrend \n');
return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = sxa_struct
  s.involid         = '';
  s.cvmstem         = '';
  s.ccvmstem        = '';
  s.maskid          = 0;
  s.thresh          = 0.5;
  s.keepmean        = 0;
  s.keeptrend       = 0;
  s.synth           = 0;
return;

%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = sxa_struct;
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

      case '-cvm',
        arg1check(flag,narg,ninputargs);
        s.cvmstem = inputargs{narg};
        narg = narg + 1;

      case '-ccvm',
        arg1check(flag,narg,ninputargs);
        s.ccvmstem = inputargs{narg};
        narg = narg + 1;

      case '-mask',
        arg1check(flag,narg,ninputargs);
        s.maskid = inputargs{narg};
        narg = narg + 1;

      case {'-thresh'} 
        arg1check(flag,narg,ninputargs);
        s.thresh = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-synth'}
        s.synth = 1;

      case {'-keepmean'}
        s.keepmean = 1;

      case {'-keeptrend'}
        s.keeptrend = 1;
  
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
%% Check argument consistency, etc %%%
function s = check_params(s)

  fprintf(1,'Checking Parameters\n');

  if(isempty(s.involid)) 
    fprintf(2,'ERROR: No input volume specified\n');
    s=[]; return;
  end

  if(isempty(s.cvmstem)) 
    fprintf(2,'ERROR: No cvm output stem specified\n');
    s=[]; return;
  end

  if(~isempty(s.ccvmstem) & isempty(s.maskid)) 
    fprintf(2,'ERROR: need mask for complementary cvm \n');
    s=[]; return;
  end





return;
%--------------------------------------------------%


