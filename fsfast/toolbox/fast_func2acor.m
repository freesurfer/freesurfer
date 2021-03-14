function r = fast_func2acor(varargin)
% r = fast_func2acor(varargin)


%
% fast_func2acor.m
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
version = 'fast_func2acor.m @FS_VERSION@';
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


print_main_struct(s,1);

s.ntrs = s.ntrs - s.nskip;

% Construct High-Pass Filter %
if(~isempty(s.hpfcutoff))
  HPF = eye(s.ntrs) - fast_mkgausmtx(s.hpfcutoff/s.tr,s.ntrs);
else
  fprintf('INFO: not using a high-pass filter\n');
  HPF = [];
end

if(~isempty(s.maskvolid))
  mask = fmri_ldbvolume(s.maskvolid);
  if(isempty(mask))
    fprintf('ERROR: could not load mask %s\n',s.maskvolid);
    return;
  end
else
  mask = [];
end

if(s.pforder >= 0)
  X = fast_polytrendmtx(1,s.ntrs,1,s.pforder);
  E = eye(s.ntrs) - X*inv(X'*X)*X';
else
  E = [];
end

minrho_mask = zeros(s.nrows*s.ncols,s.nslices);

nmasktot = 0;
yacorsum = 0;
yacorsumabs = 0;
ytwfsum = 0;
nminrhotot = 0;
for slice = s.firstslice:s.lastslice
  fprintf('%d ',slice);

  if(~isempty(mask))
    imask = find(mask(slice+1,:,:)==1);
    inotmask = find(mask(slice+1,:,:)==0);
  else
    imask = 1:s.nrows*s.ncols;
    inotmask = [];
  end

  if(isempty(imask)) continue; end
  nmask = length(imask);

  funcfile = sprintf('%s_%03d.%s',s.funcvolid,slice,s.funcbext);
  y = fmri_ldbfile(funcfile);
  if(s.nskip > 0) y = y(:,:,s.nskip+1:s.ntrs+s.nskip); end
  y = reshape(y, [s.nrows*s.ncols s.ntrs])'; %'
  if(~isempty(E))   y = E*y; end
  if(~isempty(HPF)) y = HPF*y; end

  yacor = fast_acorr(y);

  if(~isempty(s.acorvolid))
    tmp = reshape(yacor', [s.nrows s.ncols s.ntrs]); %'
    acorfile = sprintf('%s_%03d.bfloat',s.acorvolid,slice);
    fmri_svbfile(tmp,acorfile);
  end

  if(s.minrho > 0)
    indminrho = find(abs(yacor(2,:)) > s.minrho );
    minrho_mask(indminrho,slice+1) = 1;
    minrho_mask(inotmask,slice+1)  = 0;
    indminrho = find(minrho_mask(:,slice+1));
    nminrho = length(indminrho);
    fprintf('N=%d,  rho > %g\n',nminrho,s.minrho);
  else
    indminrho = imask;
    nminrho = length(indminrho);
  end
  nminrhotot = nminrhotot + nminrho;

  yacorsum = yacorsum + sum(yacor(:,indminrho),2);
  yacorsumabs = yacorsumabs + sum(abs(yacor(:,indminrho)),2);
  ytwfsum = ytwfsum + sum(y(:,imask),2);
  nmasktot = nmasktot + nmask;

end
fprintf('\n');

yacormean = yacorsum/nminrhotot;

if(~isempty(s.acormeanfile))
  fid = fopen(s.acormeanfile,'w');
  fprintf(fid,'%f\n',yacormean);
  fclose(fid);
end

if(~isempty(s.acormeanabsfile))
  yacormeanabs = yacorsumabs/nmasktot;
  fid = fopen(s.acormeanfile,'w');
  fprintf(fid,'%f\n',yacormean);
  fclose(fid);
end

if(~isempty(s.meantwffile))
  ytwfmean = ytwfsum/nmasktot;
  fid = fopen(s.acormeanfile,'w');
  fprintf(fid,'%f\n',ytwfmean);
  fclose(fid);
end

% Construct a whitening filter
c = zeros(s.ntrs,1);
c(1) = 1;
tmp = zeros(size(yacormean));
tmp(1:round(s.ntrs/2)) = yacormean(1:round(s.ntrs/2));
N = toeplitz(c,tmp);
W = inv(N);

% Construct final filter
if(~isempty(HPF)) WF = W*HPF;
else              WF = W;
end

if(~isempty(s.whtfilterfile))
  save(s.whtfilterfile,'WF');
end


% Whiten and save the output %
if(~isempty(s.whtvolid))

  fprintf('Whitening\n');
  for slice = s.firstslice:s.lastslice
    fprintf('slice %d\n',slice);

    funcfile = sprintf('%s_%03d.%s',s.funcvolid,slice,s.funcbext);
    y = fmri_ldbfile(funcfile);
    y = reshape(y, [s.nrows*s.ncols s.ntrs])'; %'
    if(~isempty(E)) y = E*y; end

    if(s.minrho > 0)
      indminrho = find(minrho_mask(:,slice+1));
    if(~isempty(indminrho))
        y(:,indminrho) = WF*y(:,indminrho);
      end
    else
      y = WF*y;
    end

    tmp = reshape(y', [s.nrows s.ncols s.ntrs]); %'
    whtfile = sprintf('%s_%03d.bfloat',s.whtvolid,slice);
    fmri_svbfile(tmp,whtfile);
  end

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
  fprintf(1,'  fast_func2acor\n');
  fprintf(1,'     -f    funcvol \n');
  fprintf(1,'     -a    acorvol \n');
  fprintf(1,'     -m    maskvol \n');
  fprintf(1,'     -am   acormean \n');
  fprintf(1,'     -ama  acormeanabs \n');
  fprintf(1,'     -mtwf mean temporal waveform\n');
  fprintf(1,'     -hpf  hpfcutoff \n');
  fprintf(1,'     -polyfit order \n');
  fprintf(1,'     -wf  whitening filter (fname.mat) \n');
  fprintf(1,'     -wht whtvol \n');
  fprintf(1,'     -minrho rho : min acor coef for whitening \n');
  fprintf(1,'     -tr  TR \n');
  fprintf(1,'     -nslices    nslices\n');
  fprintf(1,'     -firstslice s0\n');
  fprintf(1,'     -lastslice  sL\n');

return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = main_struct
  s.funcvolid      = '';
  s.funcbext       = '';
  s.acorvolid      = '';
  s.maskvolid      = '';
  s.maskbext       = '';
  s.whtvolid       = '';
  s.minrho         =  0;
  s.acormeanfile   = '';
  s.acormeanabsfile   = '';
  s.meantwffile    = '';
  s.hpfcutoff      = [];
  s.pforder        = -1;
  s.whtfilterfile  = ''; 
  s.tr             = [];
  s.nskip          = 0;
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

    case {'-a','-acorvol'}
      arg1check(flag,narg,ninputargs);
      s.acorvolid = inputargs{narg};
      narg = narg + 1;

    case {'-am','-acormean'}
      arg1check(flag,narg,ninputargs);
      s.acormeanfile = inputargs{narg};
      narg = narg + 1;

    case {'-ama','-acormeanabs'}
      arg1check(flag,narg,ninputargs);
      s.acormeanabsfile = inputargs{narg};
      narg = narg + 1;

    case {'-mtwf'}
      arg1check(flag,narg,ninputargs);
      s.meantwffile = inputargs{narg};
      narg = narg + 1;

    case {'-minrho'}
      arg1check(flag,narg,ninputargs);
      s.minrho = sscanf(inputargs{narg},'%f',1);
      narg = narg + 1;

    case {'-wht','-whtvol'}
      arg1check(flag,narg,ninputargs);
      s.whtvolid = inputargs{narg};
      narg = narg + 1;

    case {'-m','-maskvol'}
      arg1check(flag,narg,ninputargs);
      s.maskvolid = inputargs{narg};
      narg = narg + 1;

    case {'-wf'}
      arg1check(flag,narg,ninputargs);
      s.whtfilterfile = inputargs{narg};
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

    case {'-nskip'}
      arg1check(flag,narg,ninputargs);
      s.nskip = sscanf(inputargs{narg},'%d',1);
      narg = narg + 1;

    case {'-firstslice','-fs'}
      arg1check(flag,narg,ninputargs);
      s.firstslice = sscanf(inputargs{narg},'%d',1);
      narg = narg + 1;

    case {'-lastslice','-ls'}
      arg1check(flag,narg,ninputargs);
      s.lastslice = sscanf(inputargs{narg},'%d',1);
      narg = narg + 1;

    case {'-nslices','-ns'}
      arg1check(flag,narg,ninputargs);
      s.nslices = sscanf(inputargs{narg},'%d',1);
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

  %if(isempty(s.acorvolid))
  %  fprintf(2,'ERROR: No output volume specified\n');
  %  s=[]; return;
  %end

  if(isempty(s.tr))
    fprintf(2,'ERROR: No TR specified\n');
    s=[]; return;
  end

  if(~isempty(s.maskvolid))
    [nrows ncols ntp fs ns endian s.maskbext] = fmri_bfiledim(s.maskvolid);
  end

  [nrows ncols ntp fs ns endian s.funcbext] = fmri_bfiledim(s.funcvolid);
  [nslices s.nrows s.ncols s.ntrs] = fmri_bvoldim(s.funcvolid);

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

return;

%--------------------------------------------------%
%% Print data structure
function s = print_main_struct(s,fid)
  if(nargin == 1) fid = 1; end

  fprintf(fid,'FuncVol %s\n',s.funcvolid);
  fprintf(fid,'AcorVol %s\n',s.acorvolid);
  fprintf(fid,'TR  %f\n',s.tr);
  if(~isempty(s.hpfcutoff))
    fprintf(fid,'HPF Cuttoff  %f\n',s.hpfcutoff);
  else
    fprintf(fid,'No HPF used\n');
  end
  fprintf(fid,'PolyFit Order %d\n',s.pforder);
  fprintf(fid,'Min Rho       %g\n',s.minrho);
  fprintf(fid,'nskip        %d\n',s.nskip);
  fprintf(fid,'firstslice   %d\n',s.firstslice);
  fprintf(fid,'lastslice    %d\n',s.lastslice);
  fprintf(fid,'nslices      %d\n',s.nslices);

return;
%--------------------------------------------------%


