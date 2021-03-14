function r = fast_cmpdesign_nsd(varargin)
% r = fast_cmpdesign_nsd
%
% Compares two designs based on non-schedule dependent parameters.
%


%
% fast_cmpdesign_nsd.m
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

version = 'fast_cmpdesign_nsd.m @FS_VERSION@';
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

% Check that the output file can be opened %
fid = fopen(s.outfile,'w');
if(fid == -1)
  fprintf('ERROR: could not open %s\n',s.outfile);
  return;
end
fclose(fid);

if(isempty(s.seed)) s.seed = sum(100*clock); end;
if(s.seed < 0)      s.seed = sum(100*clock); end;
randn('state',s.seed); 
rand('state',s.seed); 

print_main_struct(s);

RunDuration = s.TR * s.Ntrs;

A1 = nsd2amtx(s.nsd1, s.Nc);
[F1 LPF1 HPF1 F01 XF1] = nsd2filter(s.nsd1, s.TR, s.Ntrs);
if(~isempty(F1)) F1tF1 = F1'*F1; %'
else             F1tF1 = []; 
end
Ncvm1 = [];
if(~isempty(s.nsd1.ar1rho)) Ncvm1 = fast_ar1mtx(s.nsd2.ar1rho,s.Ntrs); end

A2 = nsd2amtx(s.nsd2, s.Nc);
[F2 LPF2 HPF2 F02 XF2] = nsd2filter(s.nsd2, s.TR, s.Ntrs);
if(~isempty(F2)) F2tF2 = F2'*F2; %'
else             F2tF2 = []; 
end
Ncvm2 = [];
if(~isempty(s.nsd2.ar1rho)) Ncvm2 = fast_ar1mtx(s.nsd2.ar1rho,s.Ntrs); end

Xmntrnd1 = nsd2trendmtx(s.nsd1,s.Ntrs);
Xmntrnd2 = nsd2trendmtx(s.nsd2,s.Ntrs);

if(strcmp(s.nsd1.ermfunc,'fir'))
  TER1 = s.nsd1.ermparams(1);
  TPS1 = s.nsd1.ermparams(2);
  TW1  = s.nsd1.ermparams(3);
else
  TPS1 = 0;
  TER1 = s.nsd1.ermparams(3);
  TW1  = s.nsd1.ermparams(4);
end
if(strcmp(s.nsd2.ermfunc,'fir'))
  TER2 = s.nsd2.ermparams(1);
  TPS2 = s.nsd2.ermparams(2);
  TW2  = s.nsd2.ermparams(3);
else
  TPS2 = 0;
  TER2 = s.nsd2.ermparams(3);
  TW2  = s.nsd2.ermparams(4);
end

tic;
parlist = fast_randschedule(s.Npc,s.Tpc,RunDuration,s.tPreScan,s.Nruns,s.TR);

if(~isempty(s.svnthpar))
  par = parlist(:,:,s.svnthpar);
  fmri_svpar(par,s.parfile,[],RunDuration);
end

W1 = [];
W2 = [];

for run = 1:s.Nruns
  par = parlist(:,:,run);

  % ----------- Condition 1 ------------------- %
  if(~isempty(s.nsd1.nonlin))
    W1 = fast_wparnonlin(par,s.Tpc,s.nsd1.nonlin(1),s.nsd1.nonlin(2));
  end
  Xfir1 = fast_par2Xfir(par,s.Ntrs,s.TR,TER1,TPS1,TW1,W1);
  if(~isempty(A1)) Xtask1 = Xfir1 * A1; 
  else             Xtask1 = Xfir1;
  end
  X1 = [Xtask1 Xmntrnd1 XF1];
  if(~isempty(F1)) 
    X1tmp = F1*X1;
    tmp1 = inv(X1tmp'*X1tmp)*X1tmp'*F1;
  else             
    X1tmp = X1;
    tmp1 = inv(X1tmp'*X1tmp)*X1tmp';
  end
  %tmp1 = inv(X1tmp'*X1tmp)*X1tmp';
  if(~isempty(Ncvm1)) beta1cvm = tmp1*Ncvm1*tmp1'; %'
  else                beta1cvm = tmp1*tmp1'; %'
  end
  ind1task = 1:size(Xtask1,2);
  beta1cvmdiag = diag(beta1cvm(ind1task,ind1task));

  eff1(run,1)            = 1.0/sum(beta1cvmdiag);
  NoiseSDMult1           = sqrt(beta1cvmdiag);
  AvgNoiseSDMult1(run,1) = mean(NoiseSDMult1);
  StdNoiseSDMult1(run,1) = std(NoiseSDMult1);

  % ----------- Condition 2 ------------------- %
  if(~isempty(s.nsd2.nonlin))
    W2 = fast_wparnonlin(par,s.Tpc,s.nsd2.nonlin(1),s.nsd2.nonlin(2));
  end
  Xfir2 = fast_par2Xfir(par,s.Ntrs,s.TR,TER2,TPS2,TW2,W2);
  if(~isempty(A2)) Xtask2 = Xfir2 * A2; 
  else             Xtask2 = Xfir2;
  end
  X2 = [Xtask2 Xmntrnd2 XF2];
  if(~isempty(F2)) X2tmp = F2*X2;
  else             X2tmp = X2;
  end
  if(~isempty(F2)) 
    X2tmp = F2*X2;
    tmp2 = inv(X2tmp'*X2tmp)*X2tmp'*F2;
  else             
    X2tmp = X2;
    tmp2 = inv(X2tmp'*X2tmp)*X2tmp';
  end
  %tmp2 = inv(X2tmp'*X2tmp)*X2tmp';
  if(~isempty(Ncvm2)) beta2cvm = tmp2*Ncvm2*tmp2'; %'
  else                beta2cvm = tmp2*tmp2'; %'
  end

  ind2task = 1:size(Xtask2,2);
  beta2cvmdiag = diag(beta2cvm(ind2task,ind2task));

  eff2(run,1)            = 1.0/sum(beta2cvmdiag);
  NoiseSDMult2           = sqrt(beta2cvmdiag);
  AvgNoiseSDMult2(run,1) = mean(NoiseSDMult2);
  StdNoiseSDMult2(run,1) = std(NoiseSDMult2);

end

tmp = [eff1 eff2 AvgNoiseSDMult1 AvgNoiseSDMult2 ...
	    StdNoiseSDMult1 StdNoiseSDMult2];

fid = fopen(s.outfile,'w');
if(fid == -1)
  fprintf('ERROR: could not open %s\n',s.outfile);
  return;
end

fmt = repmat('%8.4f ',[1 size(tmp,2)]);
fmt = [fmt '\n'];
fprintf(fid,fmt,tmp'); %'
fclose(fid);

if(s.verbose) fprintf(fmt,tmp'); end %'
fprintf(fmt,mean(tmp)');%'
fprintf(fmt,max(tmp)');%'

fprintf(1,'fast_cmpdesign_nsd: Done %g\n',toc);

if(~s.rtdata) 
  r = 0;
else
  r = tmp;
end

return;
%----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%-----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%


%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = main_struct;
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

      case {'-o','-outfile'}
        arg1check(flag,narg,ninputargs);
        s.outfile = inputargs{narg};
        narg = narg + 1;

      case {'-svnthpar'}
        argNcheck(flag,narg,ninputargs,2);
        s.svnthpar = sscanf(inputargs{narg},'%d',1);
        s.parfile = inputargs{narg+1};
        narg = narg + 2;

      case '-seed',
        arg1check(flag,narg,ninputargs);
        s.seed = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-nruns',
        arg1check(flag,narg,ninputargs);
        s.Nruns = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-nc',
        arg1check(flag,narg,ninputargs);
        s.Nc = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-npc',
        if(isempty(s.Nc)) 
  	  fprintf('ERROR: need -nc before -npc\n');
          s = []; return;
        end
        argNcheck(flag,narg,ninputargs,s.Nc);
        for n = 1:s.Nc
  	  s.Npc(n) = sscanf(inputargs{narg+n-1},'%d');
        end
        narg = narg + s.Nc;

      case '-tpc',
        if(isempty(s.Nc)) 
  	  fprintf('ERROR: need -nc before -tpc\n');
          s = []; return;
        end
        argNcheck(flag,narg,ninputargs,s.Nc);
        for n = 1:s.Nc
  	  s.Tpc(n) = sscanf(inputargs{narg+n-1},'%f');
        end
        narg = narg + s.Nc;

      case '-ntrs',
        arg1check(flag,narg,ninputargs);
        s.Ntrs = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-tr',
        arg1check(flag,narg,ninputargs);
        s.TR = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-tprescan'}
        arg1check(flag,narg,ninputargs);
        s.tPreScan = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-nonlin'}
        argNcheck(flag,narg,ninputargs,2);
        s.nsd1.nonlin(1) = sscanf(inputargs{narg},'%f',1);
        s.nsd1.nonlin(2) = sscanf(inputargs{narg+1},'%f',1);
        s.nsd2.nonlin = s.nsd1.nonlin;
        narg = narg + 2;
      case {'-nonlin1'}
        argNcheck(flag,narg,ninputargs,2);
        s.nsd1.nonlin(1) = sscanf(inputargs{narg},'%f',1);
        s.nsd1.nonlin(2) = sscanf(inputargs{narg+1},'%f',1);
        narg = narg + 2;
      case {'-nonlin2'}
        argNcheck(flag,narg,ninputargs,2);
        s.nsd2.nonlin(1) = sscanf(inputargs{narg},'%f',1);
        s.nsd2.nonlin(2) = sscanf(inputargs{narg+1},'%f',1);
        narg = narg + 2;

      case {'-erm'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd1.ermfunc = inputargs{narg};
        nparams = nparams_ermfunc(s.nsd1.ermfunc);
        argNcheck(flag,narg,ninputargs,nparams+1);
        for n = 1:nparams,
          s.nsd1.ermparams(n) = sscanf(inputargs{narg+n},'%f',nparams);
        end
        s.nsd2.ermfunc   = s.nsd1.ermfunc;
        s.nsd2.ermparams = s.nsd1.ermparams;
        narg = narg + nparams + 1;
      case {'-erm1'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd1.ermfunc = inputargs{narg};
        nparams = nparams_ermfunc(s.nsd1.ermfunc);
        argNcheck(flag,narg,ninputargs,nparams+1);
        s.nsd1.ermparams = [];
        for n = 1:nparams,
          s.nsd1.ermparams(n) = sscanf(inputargs{narg+n},'%f',nparams);
        end
        narg = narg + nparams + 1;
      case {'-erm2'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd2.ermfunc = inputargs{narg};
        nparams = nparams_ermfunc(s.nsd2.ermfunc);
        argNcheck(flag,narg,ninputargs,nparams+1);
        s.nsd2.ermparams = [];
        for n = 1:nparams,
          s.nsd2.ermparams(n) = sscanf(inputargs{narg+n},'%f',nparams);
        end
        narg = narg + nparams + 1;

      case {'-lpf'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd1.LPFCutoff = sscanf(inputargs{narg},'%f',1);
        s.nsd2.LPFCutoff = s.nsd1.LPFCutoff;
        narg = narg + 1;
      case {'-lpf1'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd1.LPFCutoff = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;
      case {'-lpf2'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd2.LPFCutoff = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-hpf'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd1.HPFCutoff = sscanf(inputargs{narg},'%f',1);
        s.nsd2.HPFCutoff = s.nsd1.HPFCutoff;
        narg = narg + 1;
      case {'-hpf1'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd1.HPFCutoff = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;
      case {'-hpf2'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd2.HPFCutoff = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-fitfilter'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd1.fitfilterpve = sscanf(inputargs{narg},'%f',1);
        s.nsd2.fitfilterpve = s.nsd1.fitfilterpve;
        narg = narg + 1;
      case {'-fitfilter1'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd1.fitfilterpve = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;
      case {'-fitfilter2'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd2.fitfilterpve = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-redfilter'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd1.redfilterpve = sscanf(inputargs{narg},'%f',1);
        s.nsd2.redfilterpve = s.nsd1.redfilterpve;
        narg = narg + 1;
      case {'-redfilter1'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd1.redfilterpve = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;
      case {'-redfilter2'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd2.redfilterpve = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-ar1rho'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd1.ar1rho = sscanf(inputargs{narg},'%f',1);
        s.nsd2.ar1rho = s.nsd1.ar1rho;
        narg = narg + 1;
      case {'-ar1rho1'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd1.ar1rho = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;
      case {'-ar1rho2'}
        argNcheck(flag,narg,ninputargs,1);
        s.nsd2.ar1rho = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-meanfit'}
        s.nsd1.meanfit = 1;
        s.nsd2.meanfit = 1;
      case {'-meanfit1'}
        s.nsd1.meanfit = 1;
      case {'-meanfit2'}
        s.nsd2.meanfit = 1;

      case {'-trendfit'}
        s.nsd1.trendfit = 1;
        s.nsd2.trendfit = 1;
      case {'-trendfit1'}
        s.nsd1.trendfit = 1;
      case {'-trendfit2'}
        s.nsd2.trendfit = 1;

      case '-rtdata',
        s.rtdata = 1;

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
%% Check that there is at least N or more argument %%
function argNcheck(flag,nflag,ninputargs,N)
  if(nflag + N - 1 > ninputargs ) 
    fprintf(1,'ERROR: Flag %s needs %d arguments\n',flag,N);
    error;
  end
return;

%--------------------------------------------------%
%% Check that there is at least one more argument %%
function arg1check(flag,nflag,nmax)
  if(nflag>nmax) 
    fprintf(1,'ERROR: Flag %s needs one argument\n',flag);
    error;
  end
return;

%--------------------------------------------------%
%% Check argument consistency, etc %%%
function s = check_params(s)
  fprintf(1,'Checking Parameters\n');

  if(isempty(s.outfile))
    fprintf('ERROR: no output file specfied\n');
    s = []; return;
  end
  if(isempty(s.Nruns))
    fprintf('ERROR: -nruns not specfied\n');
    s = []; return;
  end
  if(isempty(s.Nc))
    fprintf('ERROR: -nc not specfied\n');
    s = []; return;
  end
  if(isempty(s.Npc))
    fprintf('ERROR: -npc not specfied\n');
    s = []; return;
  end
  if(isempty(s.Tpc))
    fprintf('ERROR: -tpc not specfied\n');
    s = []; return;
  end
  if(isempty(s.Ntrs))
    fprintf('ERROR: -ntrs not specfied\n');
    s = []; return;
  end
  if(isempty(s.TR))
    fprintf('ERROR: -tr not specfied\n');
    s = []; return;
  end
  s.nsd1 = check_nsd_params(s.nsd1,1);
  if(isempty(s.nsd1)) s = []; return; end
  s.nsd2 = check_nsd_params(s.nsd2,1);
  if(isempty(s.nsd2)) s = []; return; end

  if(~isempty(s.svnthpar))
    if(s.svnthpar > s.Nruns)
      fprintf('ERROR: svnthpar = %d > Nruns = %d\n',s.svnthpar,s.Nruns);
      s = []; return;
    end
  end

return;
%--------------------------------------------------%
function nsd = check_nsd_params(nsd,nsdid)

  if(isempty(nsd.ermfunc))
    fprintf('ERROR: -erm%d not specfied\n',nsdid);
    nsd = []; return;
  end
  if(~isempty(nsd.redfilterpve) & ~isempty(nsd.fitfilterpve))
    fprintf('ERROR: nsd %d, cannot -est and -fit\n',nsdid);
    nsd = []; return;
  end

return
%--------------------------------------------------%
%% Print main data structure
function print_main_struct(s)

  fprintf('outfile = %s\n',s.outfile);
  if(~isempty(s.svnthpar))
    fprintf('svnthpar = %d\n',s.svnthpar);
    fprintf('parfile  = %s\n',s.parfile);
  end
  fprintf('seed    = %g\n',s.seed);
  fprintf('Nruns   = %d\n',s.Nruns);
  fprintf('Nc      = %d\n',s.Nc);
  fprintf('Npc   = '); fprintf('%d ',s.Npc); fprintf('\n');
  fprintf('Tpc   = '); fprintf('%f ',s.Tpc); fprintf('\n');
  fprintf('Ntrs  = %d\n',s.Ntrs);
  fprintf('TR    = %d\n',s.TR);
  fprintf('tPreScan = %d\n',s.tPreScan);
  fprintf('NSD1 -------------------------\n');
  print_nsd_struct(s.nsd1)
  fprintf('NSD2 -------------------------\n');
  print_nsd_struct(s.nsd2)

return;
%--------------------------------------------------%
%% Print nsd data structure
function print_nsd_struct(nsd)

  fprintf('ERM %s ',nsd.ermfunc); fprintf('%g ',nsd.ermparams); fprintf('\n');
  if(~isempty(nsd.nonlin))
    fprintf('NonLin  %g %g\n',nsd.nonlin(1),nsd.nonlin(2));
  else
    fprintf('NonLin  None\n');
  end

  if(~isempty(nsd.LPFCutoff)) 
    fprintf('LPFCutoof  %g\n',nsd.LPFCutoff); 
  else
    fprintf('LPFCutoff  None\n');
  end

  if(~isempty(nsd.HPFCutoff)) 
    fprintf('HPFCuttof  %g\n',nsd.HPFCutoff); 
  else
    fprintf('LPFCutoff  None\n');
  end

  if(~isempty(nsd.fitfilterpve)) 
    fprintf('FitFilterPVE  %g\n',nsd.fitfilterpve); 
  else
    fprintf('FitFilterPVE  None\n');
  end

  if(~isempty(nsd.redfilterpve)) 
    fprintf('RedFilterPVE  %g\n',nsd.redfilterpve); 
  else
    fprintf('RedFilterPVE  None\n');
  end

  if(~isempty(nsd.ar1rho)) 
    fprintf('AR1 %g\n',nsd.ar1rho); 
  else 
    fprintf('AR1 None\n');
  end

  fprintf('MeanFit  %d\n',nsd.meanfit);
  fprintf('TrendFit %d\n',nsd.trendfit);

return;
%--------------------------------------------------%
%------------- Print Usage ---------------------%
function print_usage(dummy)

  fprintf(1,'USAGE:\n');
  fprintf(1,' fast_cmpdesign_nsd\n');
  fprintf(1,'     -o outfile\n');

  fprintf(1,'   \n');
  fprintf(1,'   Schedule-Dependent Parameters \n');
  fprintf(1,'   \n');
  fprintf(1,'     -nruns N : number of schedules to use\n');
  fprintf(1,'     -nc  N : number of conditions\n');
  fprintf(1,'     -npc N1 N2 ... : number of pres per cond  \n');
  fprintf(1,'     -tpc T1 T2 ... : time of pres per cond  \n');
  fprintf(1,'     -ntrs N  : number of time points\n');
  fprintf(1,'     -tr   TR : TR in seconds\n');
  fprintf(1,'     -tprescan T : duration before scanning onset\n');

  fprintf(1,'   \n');
  fprintf(1,'   Schedule-Independent Parameters \n');
  fprintf(1,'   \n');
  fprintf(1,'     -nonlin<1,2> tau attnmax : non-linear SOA effect \n');
  fprintf(1,'     -erm<1,2> func params : event response model \n');
  fprintf(1,'     -lpf<1,2> cutoff    : LPF cutoff in seconds\n');
  fprintf(1,'     -hpf<1,2> cutoff    : HPF cutoff in seconds\n');
  fprintf(1,'     -fitfilter<1,2> pve : implement filter as regressor \n');
  fprintf(1,'     -redfilter<1,2> pve : reduce dim of filter with SVD\n');
  fprintf(1,'     -ar1rho<1,2> rho    : unmodeled ar1 noise \n');
  fprintf(1,'     -meanfit<1,2>       : add regressor to fit mean \n');
  fprintf(1,'     -trendfit<1,2>      : add regressor to fit trend \n');
  fprintf(1,'   \n');
  fprintf(1,'  Possible ERMs and their parameters: \n');
  fprintf(1,'    fir       TER   tPreStim TimeWindow\n');
  fprintf(1,'    gamma     delay dispersion tres tw\n');
  fprintf(1,'    gamma+der delay dispersion tres tw\n');
  fprintf(1,'   \n');

  fprintf(1,'  Other Options \n');
  fprintf(1,'   \n');
  fprintf(1,'     -svnthpar n parfile  : save the nth par in parfile\n');

  fprintf(1,'   \n');

return
%--------------------------------------------------%
%% Main data structure
function s = main_struct
  s.seed        = [];
  s.Nruns       = [];
  s.Nc          = [];
  s.Npc         = [];
  s.Tpc         = [];
  s.Ntrs        = [];
  s.TR          = [];
  s.tPreScan    = 0;
  s.nsd1        = nsd_struct;
  s.nsd2        = nsd_struct;
  s.outfile     = [];
  s.svnthpar    = [];
  s.parfile     = [];
  s.rtdata      = 0;
  s.verbose     = 0;
return;
%--------------------------------------------------%
%% Non-schedule dependent parameters
function nsd = nsd_struct
  nsd.nonlin        = []; % [tau attnmax tres tw]
  nsd.ermfunc       = '';
  nsd.ermparams     = [];
  nsd.LPFCutoff     = [];
  nsd.HPFCutoff     = [];
  nsd.meanfit       =  0;
  nsd.trendfit      =  0;
  nsd.fitfilterpve  = [];
  nsd.redfilterpve  = [];
  nsd.ar1rho        = [];
return;
%--------------------------------------------------%
function nparams = nparams_ermfunc(ermfunc)
  nparams = -1;
  switch(ermfunc)
  case 'fir',       nparams = 3; return;
  case 'gamma',     nparams = 4; return;
  case 'gamma+der', nparams = 4; return;
  otherwise
    fprintf('ERROR: ERM %s unrecognized\n',ermfunc);
    return;
  end
return
%--------------------------------------------------%
function [F, LPF, HPF, F0, XF] = nsd2filter(nsd, TR, Ntp)
  F   = []; LPF = []; HPF = []; F0  = []; XF  = [];
  if(isempty(nsd.LPFCutoff) & isempty(nsd.HPFCutoff)) return; end

  % Low-pass Filter %
  if(~isempty(nsd.LPFCutoff))
    LPF = fast_mkgausmtx(nsd.LPFCutoff/TR,Ntp);
  else
    LPF = eye(Ntp);
  end

  % High-pass Filter %
  if(~isempty(nsd.HPFCutoff))
    HPF = eye(Ntp) - fast_mkgausmtx(nsd.HPFCutoff/TR,Ntp);
  else
    HPF = eye(Ntp);
  end

  F = LPF*HPF;
  F0 = F;

  % Reduce dimemsionality of filter based on SVD %
  if(~isempty(nsd.redfilterpve))
     [U S tmp] = svd(F);
     ds = diag(S);
     pve = 100*cumsum(ds)/sum(ds);
     ind = find(pve <= nsd.redfilterpve);
     F = U(:,ind) * S(ind,ind) * (U(:,ind)'); %'
  end

  % Get the principal components of the null space of the filter
  % for inclusion as a regressor
  if(~isempty(nsd.fitfilterpve))
     FNull = eye(size(F)) - F;
     [U S tmp] = svd(FNull);
     ds = diag(S);
     pve = 100*cumsum(ds)/sum(ds);
     ind = find(pve <= nsd.fitfilterpve);
     XF = U(:,ind);
     F = [];
  end

return;
%-------------------------------------------------------%
function A = nsd2amtx(nsd, Nc)

  A = [];
  if(~strcmp(nsd.ermfunc,'gamma') & ...
     ~strcmp(nsd.ermfunc,'gamma+der')) return; end

  delay = nsd.ermparams(1);
  dispersion = nsd.ermparams(2);
  tres  = nsd.ermparams(3);
  tw    = nsd.ermparams(4);

  a = fast_gamma(delay, dispersion, tres, 0, tw);
  a = reshape1d(a);

  d = [];
  if(strcmp(nsd.ermfunc,'gamma+der')) 
    d = fast_gammaderiv(delay, dispersion, tres, 0, tw);
    d = reshape1d(d);
  end

  b = [a d];

  A = fast_blockident(a,Nc,Nc);

return
%-------------------------------------------------------%
function Xmntrnd = nsd2trendmtx(nsd,Ntrs)

Xmntrnd = [];
if(nsd.meanfit)  
  Xmntrnd = ones(Ntrs,1); 
end
if(nsd.trendfit) 
  Xtrnd = [1:Ntrs]'; %'
  Xtrnd = Xtrnd - mean(Xtrnd);
  Xmntrnd = [Xmntrnd Xtrnd];
end

return
