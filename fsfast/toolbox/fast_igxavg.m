function r = fast_isxavg_fe(varargin)
% Name: fast_isxavg_fe
% Purpose: implements fixed-effects intersession averaging
%          for output of selxavg
% Author: Douglas Greve
% Questions or Comments: analysis-bugs@nmr.mgh.harvard.edu
% r = fast_isxavg_fe(varargin)


%
% fast_igxavg.m
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

if(s.synth < 0) s.synth = sum(100*clock); end
fprintf('SynthSeed = %10d\n',s.synth);

fprintf(1,'_______________________________________________\n');
fprintf(1,'Fixed Effects Averaging Parameters\n');
isxavg_fe_print_struct(s,1);
fprintf(1,'_______________________________________________\n');

ninvols = size(s.invols,1);
lastslice = s.firstslice + s.nslices - 1;

for slice = s.firstslice:lastslice

  % fprintf('Processing Slice %d\n',slice);
  fprintf('%2d ',slice);
  eVarSum = 0;
  DOFSum  = 0;
  hAvgSum = 0;
  hCovSum = 0;

  for session = 1:ninvols,
    fprintf('Session %d\n',session);
    InStem  = deblank(s.invols(session,:));
    InSA    = sprintf('%s_%03d.bfloat',InStem,slice);
    DatFile = sprintf('%s.dat',InStem);
    InHOffset  = sprintf('%s-offset_%03d.bfloat',InStem,slice);
    
    [hAvg eVar hd] = fast_ldsxabfile(InSA);
    if(s.synth ~= 0)
      hAvg = randn(size(hAvg));
    end

    if(s.pctsigch)
      hoffset = fmri_ldbfile(InHOffset);
      ind = find(hoffset == 0);
      hoffset(ind) = 10^10;
      eVar = eVar ./(hoffset.^2);
      hofftmp = repmat(hoffset,[1 1 size(hAvg,3)]);
      hAvg = hAvg./hofftmp;
    end

    hCov = hd.hCovMtx;
    eVarSum = eVarSum + eVar * hd.DOF;

    if(~s.weighted)
      hAvgSum = hAvgSum + hAvg * hd.DOF;
    else
      hAvgSum = hAvgSum + s.weights(session)*hAvg * hd.DOF;
    end

    hCovSum = hCovSum + inv(hCov);
    DOFSum  = DOFSum + hd.DOF;

  end % loop over sessions %

  hAvgGrp = hAvgSum/DOFSum;
  eVarGrp = eVarSum/DOFSum;
  hCovGrp = inv(hCovSum);
  hd.hCovMtx = hCovGrp;
  hd.DOF = DOFSum;
  [ySA dof] = fmri_sxa2sa(eVarGrp,hCovGrp,hAvgGrp,hd);

  OutSA = sprintf('%s_%03d.bfloat',s.avgvol,slice); 
  fmri_svbfile(ySA, OutSA); 

  OutDat = sprintf('%s.dat',s.avgvol);
  fmri_svdat2(OutDat,hd);

end % Loop over slices %
fprintf('\n');

fprintf(1,'fast_isavg_fe completed SUCCESSUFLLY\n\n');

return;

%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = isxavg_fe_struct;
  inputargs = varargin{1};
  ninputargs = length(inputargs);

  narg = 1;
  while(narg <= ninputargs)

    flag = deblank(inputargs{narg});
    narg = narg + 1;
    % fprintf(1,'Argument: %s\n',flag);
    if(~isstr(flag))
      flag
      fprintf(1,'ERROR: All Arguments must be a string\n');
      s = []; return;
    end

    switch(flag)

      case '-weighted',
        s.weighted = 1;

      case '-pctsigch',
        s.pctsigch = 1;

      case '-i',
        if(~s.weighted) 
          if(arg1check(flag,narg,ninputargs)) s = []; return; end
        else
          if(arg2check(flag,narg,ninputargs)) s = []; return; end
        end
        s.invols = strvcat(s.invols,inputargs{narg});
        narg = narg + 1;

        if(s.weighted) 
          n = size(s.invols,1);
          s.weights(n) = sscanf(inputargs{narg},'%f',1);
          narg = narg + 1;
        end

      case {'-firstslice', '-fs'}
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.firstslice = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-nslices', '-ns'}
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.nslices = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-avg','-o'},
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.avgvol = inputargs{narg};
        narg = narg + 1;

      case '-logfile',
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.logfile = inputargs{narg};
        narg = narg + 1;

      case '-synth',
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.synth  = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case '-nolog',
        s.nolog = 1;
        s.logfid = 0;

      case {'-monly', '-umask'}, % ignore
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        narg = narg + 1;

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
function r = arg1check(flag,nflag,nmax)
  r = 0;
  if(nflag>nmax) 
    fprintf(1,'ERROR: Flag %s needs one argument',flag);
    r = 1;
  end
return;

%--------------------------------------------------%
%% Check that there are at least two more arguments %%
function r = arg2check(flag,nflag,nmax)
  r = 0;
  if(nflag+1 > nmax) 
    fprintf(1,'ERROR: Flag %s needs two arguments',flag);
    r = 1;
  end
return;

%--------------------------------------------------%
%% Check argument consistency, etc %%%
function s = check_params(s)

  fprintf(1,'Checking Parameters\n');
  if(size(s.invols,1) < 2) 
    fprintf(1,'ERROR: must have at least 2 inputs\n');
    s = []; return;
  end

return;

%--------------------------------------------------%
%% Print Usage 
function print_usage
  fprintf(1,'USAGE:\n');
  fprintf(1,'  fast_isxavg_fe\n');
  fprintf(1,'     -i invol1 <weight> -i invol2 <weight>... \n');
  fprintf(1,'     -weighted\n');
  fprintf(1,'     -avg     avgvol \n');
  fprintf(1,'     -pctsigch \n');
  fprintf(1,'     -firstslice sliceno  \n');
  fprintf(1,'     -nslices    nslices  \n');
return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = isxavg_fe_struct
  s.invols    = '';
  s.weighted  = 0;
  s.weights   = [];
  s.firstslice = 0;
  s.nslices    = -1;
  s.avgvol = '';
  s.pctsigch = 0;
  s.logfile = '';
  s.logfid  = 1;
  s.nolog   = 0;
  s.synth   = 0;
return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Print data structure
function s = isxavg_fe_print_struct(s,fid)
  if(nargin == 1) fid = 1; end

  ninvols = size(s.invols,1);
  fprintf(fid,'ninvols    %d\n',ninvols);
  for n = 1:ninvols
    if(~s.weighted)
      fprintf(fid,'  invol      %s\n',s.invols(n,:));
    else
      fprintf(fid,'  invol      %s  %g\n',s.invols(n,:),s.weights(n));
    end
  end
  fprintf(fid,'outvol      %s\n',s.avgvol);
  fprintf(fid,'pctsigch    %d\n',s.pctsigch);
  fprintf(fid,'firstslice  %d\n',s.firstslice);
  fprintf(fid,'nslices     %d\n',s.nslices);
  fprintf(fid,'logfile     %s\n',s.logfile);
  fprintf(fid,'synth       %s\n',s.synth);

return;
%--------------------------------------------------%

