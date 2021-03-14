function r = afni2bvol(varargin)
% r = afni2bvol(varargin)


%
% afni2bvol.m
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

brik =  ReadBRIK2(s.inbrik);
if(isempty(brik)) return; end

% Permute %
brik = permute(brik,[3 2 1 4]);

fmri_svbvolume(brik,s.outvolid,size(brik),s.oext);

fprintf('afni2bvol: done\n');

return;
%----------------------------------------------------------%
%----------------------------------------------------------%
%----------------------------------------------------------%
%----------------------------------------------------------%


%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  %fprintf(1,'Parsing Arguments \n');
  s = inorm_struct;
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

     case {'-i','-inbrik'},
        arg1check(flag,narg,ninputargs);
        s.inbrik = inputargs{narg};
        narg = narg + 1;

     case {'-o','-outstem'},
        arg1check(flag,narg,ninputargs);
        s.outvolid = inputargs{narg};
        narg = narg + 1;

     case {'-oext'}
        arg1check(flag,narg,ninputargs);
        s.oext = inputargs{narg};
        narg = narg + 1;
        if(~strcmp(s.oext,'bshort') & ~strcmp(s.oext,'bfloat'))
	   fprintf('ERROR: output extension = %s, must be either\n');
	   fprintf('       bshort or bfloat\n');
           s = []; return;
        end

      case '-verbose',
        s.verbose = 1;

      case {'-monly','-umask'}, % ignore
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
%% Print Usage 
function print_usage(dummy)

  fprintf(1,'USAGE:\n');
  fprintf(1,'  afni2bvol\n');
  fprintf(1,'     -i input brik (eg, 3d+orig.BRIK)\n');
  fprintf(1,'     -o outstem \n');
  fprintf(1,'     -oext output extension (<bfloat>,bshort) \n');
return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = inorm_struct
  s.inbrik       = '';
  s.outvolid     = '';
  s.outvolext    = 'bfloat';
return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Check argument consistency, etc %%%
function s = check_params(s)

  if(isempty(s.inbrik))
    fprintf('ERROR: must specify input brik');
    s = []; return;
  end

  if(isempty(s.outvolid))
    fprintf('ERROR: must specify output stem');
    s = []; return;
  end

return;
