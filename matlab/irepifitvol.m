function r = irepifitvol(varargin)

r = 1;

cmdargs.involfile = '';
cmdargs.maskfile = '';
cmdargs.configfile = '';
cmdargs.outdir = '';
cmdargs.skip = 7;
cmdargs.nacqexclude = 0;
cmdargs.nminexclude = 0;
cmdargs.ndummies = 10;
cmdargs.ROFlip = 65; % Readout flip angle degrees
cmdargs.TBI = 2.092*1000; % ms, time bet inversions, like the TR
cmdargs.PreInv = 14; % ms, subtract this from TBS of final slice
cmdargs.TI1 = 24; % ms, time to first acq after inversion
cmdargs.InvDur = 12; %  ms, duration of inversion pulse, "Prefill"
cmdargs.SaveYHat = 0;
cmdargs.Slice1PreInv = 0; % Set to 1 to model slice 1 acqed before inv
cmdargs.T1Range = [500 50 8000]; % min delta max

%% Print useage if there are no arguments %%
if(nargin == 0)
  print_usage(cmdargs)
  return;
end

%% Parse the arguments %%
cmdargs = parse_args(cmdargs,varargin);
if(isempty(cmdargs)) return; end
cmdargs = check_params(cmdargs);
if(isempty(cmdargs)) return; end

tic;

dce = MRIread(cmdargs.involfile);
if(isempty(dce)) return; end

UseMask = 0;
if(~isempty(cmdargs.maskfile))
  mask = MRIread(cmdargs.maskfile);
  if(isempty(mask)) return; end
  UseMask = 1;
end

mkdirp(cmdargs.outdir);

nslices = dce.volsize(3);
ntp = dce.nframes;

skip = cmdargs.skip;
ndummies = cmdargs.ndummies;
ROFlipDeg = cmdargs.ROFlip; 
TBI = cmdargs.TBI; 
PreInv = cmdargs.PreInv;
InvDur = cmdargs.InvDur;
TI1 = cmdargs.TI1;

s0 = irepistructure(TBI,nslices,ndummies,ntp,skip,PreInv,InvDur,TI1,ROFlipDeg);
s0.Slice1PreInv = cmdargs.Slice1PreInv;
s0 = irepitiming(s0);
s0.nminexclude = cmdargs.nminexclude;
s0.nexclude = cmdargs.nacqexclude;
s0.T1 = [cmdargs.T1Range(1):cmdargs.T1Range(2):cmdargs.T1Range(3)];
nT1 = length(s0.T1);

fprintf('Skip = %d\n',cmdargs.skip);
fprintf('ndummies = %d\n',cmdargs.ndummies);
fprintf('ROFlip = %f\n',cmdargs.ROFlip);
fprintf('TBI = %f\n',cmdargs.TBI);
fprintf('PreInv = %f\n',cmdargs.PreInv);
fprintf('InvDur = %f\n',cmdargs.InvDur);
fprintf('TI1 = %f\n',cmdargs.TI1);
fprintf('NAcqEx = %d\n',cmdargs.nacqexclude);
fprintf('NMinEx = %d\n',cmdargs.nminexclude);
fprintf('Slice1PreInv = %d\n',cmdargs.Slice1PreInv);
fprintf('T1Range %g %g %g\n',cmdargs.T1Range);

fname = sprintf('%s/info.mat',cmdargs.outdir);
save(fname,'cmdargs','s0');

t1map = dce;
t1map.vol = zeros(dce.volsize);
rstdmap = dce;
rstdmap.vol = zeros(dce.volsize);
M0map = dce;
M0map.vol = zeros(dce.volsize);
if(cmdargs.SaveYHat)
  yhatmap = dce;
  yhatmap.vol = zeros([dce.volsize ntp]);
end

for sliceno = 1:nslices
  fprintf('sliceno = %d/%d %g\n',sliceno,nslices,toc);  
  s2 = s0;
  s2.sliceno = sliceno;
  s3 = s0;
  s3.sliceno = sliceno;
  for r = 1:size(dce.vol,1)
    for c = 1:size(dce.vol,2)
      if(UseMask & mask.vol(r,c,sliceno) < 0.5) continue; end

      y = squeeze(dce.vol(r,c,sliceno,:));
      s2.y = y;
      s2 = irepifit(s2);
      [mm imin] = min(s2.rstd);

      s3.y = y;
      s3.T1 = s2.T1(imin)+[-25:1:25];
      s3 = irepifit(s3);
      [mm imin] = min(s3.rstd);

      t1map.vol(r,c,sliceno) = s3.T1(imin);
      M0map.vol(r,c,sliceno) = s3.M0(imin);
      rstdmap.vol(r,c,sliceno) = s3.rstd(imin);

      if(cmdargs.SaveYHat)
	yhatmap.vol(r,c,sliceno,:) = s3.yhat0(:,imin);
      end

    end
  end
  fname = sprintf('%s/t1.nii.gz',cmdargs.outdir);
  MRIwrite(t1map,fname);
end

fname = sprintf('%s/t1.nii.gz',cmdargs.outdir);
MRIwrite(t1map,fname);

fname = sprintf('%s/M0.nii.gz',cmdargs.outdir);
MRIwrite(M0map,fname);

fname = sprintf('%s/rstd.nii.gz',cmdargs.outdir);
MRIwrite(rstdmap,fname);

if(cmdargs.SaveYHat)
  fname = sprintf('%s/yhat.nii.gz',cmdargs.outdir);
  MRIwrite(yhatmap,fname);
end

	
fprintf('irepifitvol done t=%g\n',toc);

return;

%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function cmdargs = parse_args(cmdargs,varargin)

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
    
   case '--i',
    arg1check(flag,narg,ninputargs);
    cmdargs.involfile = inputargs{narg};
    narg = narg + 1;
    
   case '--m',
    arg1check(flag,narg,ninputargs);
    cmdargs.maskfile = inputargs{narg};
    narg = narg + 1;
    
   case '--o',
    arg1check(flag,narg,ninputargs);
    cmdargs.outdir = inputargs{narg};
    narg = narg + 1;
    
   case '--skip',
    arg1check(flag,narg,ninputargs);
    cmdargs.skip = sscanf(inputargs{narg},'%d');
    narg = narg + 1;
    
   case '--nacqex',
    arg1check(flag,narg,ninputargs);
    cmdargs.nacqexclude = sscanf(inputargs{narg},'%d');
    narg = narg + 1;
    
   case '--nminex',
    arg1check(flag,narg,ninputargs);
    cmdargs.nminexclude = sscanf(inputargs{narg},'%d');
    narg = narg + 1;
    
   case '--ndummies',
    arg1check(flag,narg,ninputargs);
    cmdargs.ndummies = sscanf(inputargs{narg},'%d');
    narg = narg + 1;
    
   case '--roflip',
    arg1check(flag,narg,ninputargs);
    cmdargs.ROFlip = sscanf(inputargs{narg},'%f');
    narg = narg + 1;
    
   case '--t1-range',
    arg3check(flag,narg,ninputargs);
    cmdargs.T1Range(1) = sscanf(inputargs{narg},'%f');
    cmdargs.T1Range(2) = sscanf(inputargs{narg+1},'%f');
    cmdargs.T1Range(3) = sscanf(inputargs{narg+2},'%f');
    narg = narg + 3;
   
   case '--tbi',
    arg1check(flag,narg,ninputargs);
    cmdargs.TBI = sscanf(inputargs{narg},'%f');
    narg = narg + 1;
    
   case '--preinv',
    arg1check(flag,narg,ninputargs);
    cmdargs.PreInv = sscanf(inputargs{narg},'%f');
    narg = narg + 1;
    
   case '--invdur',
    arg1check(flag,narg,ninputargs);
    cmdargs.InvDur = sscanf(inputargs{narg},'%f');
    narg = narg + 1;

   case '--ti1',
    arg1check(flag,narg,ninputargs);
    cmdargs.TI1 = sscanf(inputargs{narg},'%f');
    narg = narg + 1;
    
   case '--slice1preinv',
    cmdargs.Slice1PreInv = 1; % Set to 1 to model slice 1 acqed before inv    
    
   case '--save-yhat',
    cmdargs.SaveYHat = 1;    

   case '--c',
    arg1check(flag,narg,ninputargs);
    cmdargs.configfile = inputargs{narg};
    narg = narg + 1;
    
   case '-debug',
    cmdargs.debug = 1;
    
   otherwise
    fprintf(2,'ERROR: Flag %s unrecognized\n',flag);
    s = [];
    return;
    
  end % --- switch(flag) ----- %

end % while(narg <= ninputargs)

return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Check argument consistency, etc %%%
function cmdargs = check_params(cmdargs)
  if(isempty(cmdargs.involfile))
    fprintf('ERROR: must specify an input with --i\n');
    error;
  end
  if(0 & isempty(cmdargs.configfile))
    fprintf('ERROR: must specify a config file with --c\n');
    error;
  end
  if(isempty(cmdargs.outdir))
    fprintf('ERROR: must specify an output dir with --o\n');
    error;
  end
return;

%--------------------------------------------------%
%% Check that there is at least one more argument %%
function arg1check(flag,nflag,nmax)
  if(nflag>nmax) 
    fprintf(1,'ERROR: Flag %s needs one argument',flag);
    error;
  end
return;
%--------------------------------------------------%
%% Check that there are at least two more arguments %%
function arg2check(flag,nflag,nmax)
  if(nflag > nmax-1 ) 
    fprintf(1,'ERROR: Flag %s needs two arguments',flag);
    error;
  end
return;
%--------------------------------------------------%
%% Check that there are at least two more arguments %%
function arg3check(flag,nflag,nmax)
  if(nflag > nmax-2 ) 
    fprintf(1,'ERROR: Flag %s needs three arguments',flag);
    error;
  end
return;

%------------- Print Usage ---------------------%
function print_usage(cmdargs)
  fprintf('USAGE:\n');
  fprintf('irepifitvol\n');
  fprintf(' --o outdir : output folder\n');
  fprintf(' Inputs \n');
  fprintf(' --i input : dce volume (eg, dce.nii.gz)\n');
  fprintf(' --m mask (optional, speeds processing)\n');
  fprintf(' Pulse sequence parameters \n');
  fprintf(' --skip skip : slice permutation skip factor \n');
  fprintf(' --ndummies ndummies : number of dummy scans\n');
  fprintf(' --roflip FlipAngle : Readout flip angle (degrees)\n');
  fprintf(' --tbi TBI : time between inversions, like TR (ms)\n');
  fprintf(' --preinv PreInv : subtract from TBS of last slice (ms)\n');
  fprintf(' --invdur InvDur : duration of inv pulse (ms)\n');
  fprintf(' --ti1 TI1 : time from end of inv pulse to first readout (ms)\n');
  fprintf(' --slice1preinv : assume slice1 before inversion\n');
  fprintf(' --t1-range t1min deltat1 t1max\n');
  fprintf(' Analysis parameters \n');
  fprintf(' --nacqex nacqex :  skip the first nacqex after inversion\n');
  fprintf(' --nminex nminex :  skip the nminex smallest time points\n');
  %fprintf('--c configfile ... \n');
return
%--------------------------------------------------%
