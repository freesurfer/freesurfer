function r = fast_functwf(varargin)
% r = fast_functwf(varargin)
%
% Temporal analysis of functional data.


%
% fast_functwf.m
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

tic;
version = 'fast_functwf.m @FS_VERSION@';
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

tic;

print_main_struct(s,1);

% ------------ Construct High-Pass Filter ------------------%
if(~isempty(s.hpfcutoff))
  HPF = eye(s.ntrs) - fast_mkgausmtx(s.hpfcutoff/s.tr,s.ntrs);
else
  fprintf('INFO: not using a high-pass filter\n');
  HPF = [];
  s.hpfcutoff = -1;
end

%--------------- Remove Mean and Linear Trend --------------%
if(s.pforder >= 0)
  X = fast_polytrendmtx(1,s.ntrs,1,s.pforder);
  E = eye(s.ntrs) - X*inv(X'*X)*X';
else
  E = [];
end

%--------------- Mask -----------------------------------%
if(~isempty(s.maskvolid))
  fprintf('INFO: Loading mask volume %g\n',toc);
  mask = fmri_ldbvolume(s.maskvolid);
  if(isempty(mask))
    fprintf('ERROR: could not load mask %s\n',s.maskvolid);
    r = 1; return;
  end
  imask = find(mask);
  if(isempty(imask))
    fprintf('ERROR: mask is empty\n');
    r = 1; return;
  end
  nmask = length(imask);
  fprintf('INFO: found %d voxels in mask\n',nmask);
else
  imask = [];
  nmask = -1;
end

%--------------- Input Volume --------------------------%
fprintf('INFO: Loading input volume %g\n',toc);
y = fmri_ldbvolume(s.funcvolid);
if(isempty(y))
  fprintf('ERROR: Loading volume\n');
  r = 1; return;
end
[nslices nrows ncols ntrs] = size(y);
nvoxs = prod([nslices nrows ncols]);
if(s.cutends) y([1 nslices],:,:,:) = 0; end

y = fast_voltomat2d(y);

if(~isempty(imask)) 
  if(max(imask) > nvoxs)
    fprintf('ERROR: mask has more voxels than input volume\n');
    r = 1; return;
  end
  y = y(:,imask); 
end
if(~isempty(E))   y = E*y; end
if(~isempty(HPF)) y = HPF*y; end
ny = size(y,2);

fprintf('TWF Mean and Std  %g\n',toc);
twf_mean = mean(y,2);
twf_std  = std(y,[],2);
twf_stderr  = twf_std/sqrt(ny);

fprintf('Abs TWF Mean and Std  %g\n',toc);
abstwf_mean = mean(abs(y),2);
abstwf_std  = std(abs(y),[],2);
abstwf_stderr = abstwf_std/sqrt(ny);

fprintf('Autocorrelation  %g\n',toc);
acor = fast_acorr(y);

fprintf('Autocorrelation Mean and Std  %g\n',toc);
acor_mean = mean(acor,2);
acor_std  = std(acor,[],2);
acor_stderr  = acor_std/sqrt(ny);

fprintf('Abs Autocorrelation Mean and Std  %g\n',toc);
absacor_mean = mean(abs(acor),2);
absacor_std  = std(abs(acor),[],2);
absacor_stderr = absacor_std/sqrt(ny);

clear acor;

fprintf('FFT and Power Spectrum %g\n',toc);
yfft = fft(y)/s.ntrs;
pwrspect = abs(yfft).^2;
clear yfft;

fprintf('Mean and Std Power Spectrum  %g\n',toc);
pwrspect_mean = mean(pwrspect,2);
pwrspect_std  = std(pwrspect,[],2);
pwrspect_stderr = pwrspect_std/sqrt(ny);
clear pwrspect;

% Frequency Axis, etc
freqmax = (1/s.tr)/2;  % Nyquist
deltafreq = freqmax/(s.ntrs/2);  % Measured from 0 to Nyquist
freqaxis = deltafreq * [0:s.ntrs-1]'; %'
freq_ind_zero = [round(s.ntrs/2):s.ntrs];
pwrspect_mean(freq_ind_zero) = 0; 
pwrspect_std(freq_ind_zero) = 0; 

fprintf('Temporal Covariance Matrix  %g\n',toc);
tcvm = y*y'; %'
[Tev EigVals blah] = svd(tcvm);

eigspect = diag(EigVals);
pve = 100*eigspect/sum(eigspect);
cpve = cumsum(pve);

clear y;

fprintf('Final Output  %g\n',toc);
n0 = [0:s.ntrs-1]'; %'
n1 = n0 + 1;
t  = s.tr*n0;

datamtx = zeros(s.ntrs, 1 + 49 + s.ntrs );

% The first and second columns are reserved to put 
% info about how the conditions underwhich the 
% results were computed.
datamtx(1,1) = s.pforder;
datamtx(2,1) = s.hpfcutoff;
datamtx(3,1) = ny;
datamtx(4,1) = nmask;
datamtx(5,1) = deltafreq;
datamtx(6,1) = freqmax;
datamtx(7,1) = s.cutends;
% Col 2 is reserved

% The next 98 are for the data
datamtx(:, 3) = n0;
datamtx(:, 4) = n1;
datamtx(:, 5) = t;
datamtx(:, 6) = freqaxis;

datamtx(:, 7) = twf_mean;
datamtx(:, 8) = twf_std;
datamtx(:, 9) = twf_stderr;

datamtx(:,10) = abstwf_mean;
datamtx(:,11) = abstwf_std;
datamtx(:,12) = abstwf_stderr;

datamtx(:,13) = acor_mean;
datamtx(:,14) = acor_std;
datamtx(:,15) = acor_stderr;

datamtx(:,16) = absacor_mean;
datamtx(:,17) = absacor_std;
datamtx(:,18) = absacor_stderr;

datamtx(:,19) = pwrspect_mean;
datamtx(:,20) = pwrspect_std;
datamtx(:,21) = pwrspect_stderr;

datamtx(:,22) = eigspect;
datamtx(:,23) = zeros(size(eigspect)); % reserved for expect eigspect
datamtx(:,24) = pve;
datamtx(:,25) = zeros(size(pve)); % reserved for expected pve
datamtx(:,26) = cpve;
datamtx(:,27) = zeros(size(cpve)); % reserved for expected cpve

% 28-100 - reserved

% Now put the temporal eigenvectors
datamtx(:,101:101+s.ntrs-1) = Tev;

% Write the data to a file %
ndatacols = size(datamtx,2);
fmt = repmat('%g ',[1 ndatacols]);
fmt = [fmt '\n'];

fid = fopen(s.outfile,'w');
fprintf(fid,fmt,datamtx'); %'
fclose(fid);

if(~isempty(s.tevstem))
  fprintf('Saving temporal eigenvectors to %s\n',s.tevstem);
  tmp = reshape(Tev',[1 1  size(Tev')]);
  fmri_svbvolume(tmp,s.tevstem);
end


fprintf(1,'Done %g\n',toc);

r = 0;

return;
%---\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%
%-----\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\-----%


%------------- Print Usage ---------------------%
function print_usage(dummy)

  fprintf(1,'USAGE:\n');
  fprintf(1,'  fast_functwf\n');
  fprintf(1,'     -f    funcvol \n');
  fprintf(1,'     -m    maskvol \n');
  fprintf(1,'     -hpf  hpfcutoff \n');
  fprintf(1,'     -meanfit \n');
  fprintf(1,'     -trendfit \n');
  fprintf(1,'     -o    output (default is funcvol.twf) \n');
  fprintf(1,'     -tev stem : save temp eigvects in stem_000.bfloat \n');
  fprintf(1,'     -tr  TR \n');
  fprintf(1,'\n');
  %fprintf(1,'     -nslices    nslices\n');
  %fprintf(1,'     -firstslice s0\n');
  %fprintf(1,'     -lastslice  sL\n');
  %fprintf(1,'\n');

return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = main_struct
  s.funcvolid      = '';
  s.outfile        = '';
  s.tevstem       = '';
  s.maskvolid      = '';
  s.maskbext       = '';
  s.hpfcutoff      = [];
  s.pforder        = -1;
  s.tr             = [];
  s.cutends        = 0;
  s.firstslice     = [];
  s.lastslice      = [];
  s.nslices        = [];
  s.ncols          = [];
  s.nrows          = [];
  s.ntrs           = [];
  s.debug = 0;
return;

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

    case {'-f','-funcvol'}
      arg1check(flag,narg,ninputargs);
      s.funcvolid = inputargs{narg};
      narg = narg + 1;

    case {'-m','-maskvol'}
      arg1check(flag,narg,ninputargs);
      s.maskvolid = inputargs{narg};
      narg = narg + 1;

    case {'-o'}
      arg1check(flag,narg,ninputargs);
      s.outfile = inputargs{narg};
      narg = narg + 1;

    case {'-tev'}
      arg1check(flag,narg,ninputargs);
      s.tevstem = inputargs{narg};
      narg = narg + 1;

    case {'-hpf'}
      arg1check(flag,narg,ninputargs);
      s.hpfcutoff = sscanf(inputargs{narg},'%f',1);
      narg = narg + 1;

    case {'-polyfit'}
      arg1check(flag,narg,ninputargs);
      s.pforder = sscanf(inputargs{narg},'%d',1);
      narg = narg + 1;

    case {'-tr','-TR'}
      arg1check(flag,narg,ninputargs);
      s.tr = sscanf(inputargs{narg},'%f',1);
      narg = narg + 1;

    case {'-firstslice','-fs'}
      arg1check(flag,narg,ninputargs);
      s.firstslice = sscanf(inputargs{narg},'%f',1);
      narg = narg + 1;

    case {'-cutends'}
      s.cutends = 1;

    case {'-lastslice','-ls'}
      arg1check(flag,narg,ninputargs);
      s.lastslice = sscanf(inputargs{narg},'%f',1);
      narg = narg + 1;

    case {'-nslices','-ns'}
      arg1check(flag,narg,ninputargs);
      s.nslices = sscanf(inputargs{narg},'%f',1);
      narg = narg + 1;

    % ignore these guys %
    case {'-monly','-umask'},
      arg1check(flag,narg,ninputargs);
      narg = narg + 1;

    case {'-debug','-echo'}, % ignore
      s.debug = 1;

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
%% Check that there are at least two more arguments %%
function arg2check(flag,nflag,nmax)
  if(nflag > nmax-1 ) 
    fprintf(1,'ERROR: Flag %s needs two arguments',flag);
    error;
  end
return;
%--------------------------------------------------%
%% Check argument consistency, etc %%%
function s = check_params(s)

  fprintf(1,'Checking Parameters\n');

  if(isempty(s.funcvolid))
    fprintf(2,'ERROR: No input volume specified\n');
    s=[]; return;
  end

  if(isempty(s.tr))
    fprintf(2,'ERROR: No TR specified\n');
    s=[]; return;
  end

  if(~isempty(s.maskvolid))
    [nrows ncols ntp fs ns endian s.maskbext] = fmri_bfiledim(s.maskvolid);
  end

  [nrows ncols ntp fs ns endian s.funcbext] = fmri_bfiledim(s.funcvolid);
  [nslices s.nrows s.ncols s.ntrs] = fmri_bvoldim(s.funcvolid);

  % This trap is needed for the data matrix %
  if(s.ntrs < 7)
    fprintf(2,'ERROR: ntrs = %d, at least 7 are need \n');
    s=[]; return;
  end

  if(isempty(s.firstslice)) s.firstslice = 0; end

  if(~isempty(s.lastslice) & ~isempty(s.nslices) )
    fprintf(2,'ERROR: cannot spec last slice and nslices \n');
    s=[]; return;
  end

  if(isempty(s.lastslice) & ~isempty(s.nslices))
    s.lastslice = s.nslices + s.firstslice - 1;
  end

  if(~isempty(s.lastslice) & isempty(s.nslices))
    s.nslices = s.lastslice - s.firstslice + 1;
  end

  if(isempty(s.lastslice) & isempty(s.nslices))
    s.nslices = nslices;
    s.lastslice = s.nslices + s.firstslice - 1;
  end

  if( (s.pforder > 0) & ~isempty(s.hpfcutoff))
    fprintf(2,'ERROR: cannot use polyfit and highpass filter \n');
    s=[]; return;
  end

  if(isempty(s.outfile))
    s.outfile = sprintf('%s.twf',s.funcvolid);
  end

return;

%--------------------------------------------------%
%% Print data structure
function s = print_main_struct(s,fid)
  if(nargin == 1) fid = 1; end

  fprintf(fid,'Output File %s\n',s.outfile);
  fprintf(fid,'FuncVol %s\n',s.funcvolid);
  fprintf(fid,'TR  %f\n',s.tr);

  if(~isempty(s.maskvolid))
    fprintf(fid,'Mask Volume %s\n',s.maskvolid);
  else
    fprintf(fid,'No Mask specified\n');
  end

  if(~isempty(s.hpfcutoff))
    fprintf(fid,'HPF Cuttoff  %f\n',s.hpfcutoff);
  else
    fprintf(fid,'No HPF used\n');
  end

  fprintf(fid,'PolyFit Order %d\n',s.pforder);
  fprintf(fid,'cutends       %d\n',s.cutends);
  fprintf(fid,'firstslice    %d\n',s.firstslice);
  fprintf(fid,'lastslice     %d\n',s.lastslice);
  fprintf(fid,'nslices       %d\n',s.nslices);

return;
%--------------------------------------------------%


