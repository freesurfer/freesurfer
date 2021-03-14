function r = fast_inorm(varargin)
% r = fast_inorm(varargin)
% Intensity normalizes and does some simple data hygene analysis.


%
% fast_inorm.m
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

version = 'fast_inorm.m @FS_VERSION@';
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

print_struct(s);

tic;

% --- Compute the slice-by-slice mean --- %
fprintf(1,'First Pass: Computing global mean\n');
slicemean = zeros(s.nslices,1);
nthslice = 1;
for slice = s.firstslice:s.lastslice
  fprintf(1,'%2d ',slice);
  if(rem(slice+1,10)==0) fprintf(1,'\n'); end
  fname = sprintf('%s_%03d.%s',s.invol,slice,s.involext);
  y = fmri_ldbfile(fname);
  slicemean(nthslice) = mean(reshape1d(y));
  nthslice = nthslice + 1;
end
fprintf(1,'\n'); 

% --- Compute the global mean --- %
globalmean = mean(slicemean);

% Compute the segmentation threshold %
if(isempty(s.absthresh))
  % Absolute threshold not specified, compute it %
  s.absthresh = s.relthresh * globalmean;
else
  % Relattive threshold not specified, compute it %
  % This is only for information.
  s.relthresh = globalmean/s.absthresh
end

s.absthresh_under = s.relthresh_under * globalmean;

% --- Compute the slice-by-slice average temporal waveform --- %
fprintf(1,'Second Pass: Segmenting\n');
nthslice = 1;
for slice = s.firstslice:s.lastslice
  fprintf(1,'%2d ',slice);
  if(rem(slice+1,10)==0) fprintf(1,'\n'); end

  fname = sprintf('%s_%03d.%s',s.invol,slice,s.involext);
  y = fmri_ldbfile(fname);
  y2 = reshape(y,[s.nrows*s.ncols s.ntp]);

  % Compute the temporal average at each voxel %
  mnslice = mean(y,3);

  % Segment the super-theshold voxels from the mean volume%
  iover  = find(mnslice >  s.absthresh);
  nover(nthslice,1) = length(iover);

  % Compute the spatial average at each time point %
  if(nover(nthslice,1) > 1)
    twf_over(nthslice,:) = mean(y2(iover,:));
  elseif(nover(nthslice,1) == 1)
    twf_over(nthslice,:) = y2(iover,:);
  else
    twf_over(nthslice,:) = zeros(1,s.ntp);
  end

  % Segment the sub-theshold voxels from the mean volume%
  iunder  = find(mnslice <=  s.absthresh_under);
  nunder(nthslice,1) = length(iunder);

  % Compute the spatial average at each time point %
  if(nunder(nthslice,1) > 1)
    twf_under(nthslice, :) = mean(y2(iunder,:));
  elseif (nunder(nthslice,1) == 1)
    twf_under(nthslice, :) = y2(iunder,:);
  else
    fprintf('INFO: no voxels found under threshold\n');
    twf_under(nthslice,:) = zeros(1,s.ntp);
  end

  nthslice = nthslice + 1;
end
fprintf(1,'\n'); 

%--- Deconvolution matrix for computing trends ---%
X = [ones(s.ntp,1) [0:s.ntp-1]']; %'

%-------------  SupraThreshold: -------------------------%
% --------Compute slice-based stats ------------- %
twf_over_min  = min(twf_over')';
twf_over_max  = max(twf_over')';
twf_over_mean = mean(twf_over')';
twf_over_stddev = std(twf_over')';
twf_over_absdev = abs(twf_over - repmat(twf_over_mean,[1 s.ntp]));
twf_over_meanabsdev = mean(twf_over_absdev')';
twf_over_stddev_sf = twf_over_meanabsdev;
tmp = twf_over_stddev; ind0 = find(tmp == 0); tmp(ind0) = 10^10; %Div0
twf_over_z      = twf_over_absdev./(repmat(tmp,[1 s.ntp]));
twf_over_z_max  = max(twf_over_z')';
tmp = twf_over_meanabsdev; ind0 = find(tmp == 0); tmp(ind0) = 10^10; %Div0
twf_over_snr_sf = twf_over_mean./tmp;
twf_over_z_sf   = twf_over_absdev./(repmat(tmp,[1 s.ntp]));
nout_over = zeros(s.nslices,1);
[i j] = find(twf_over_z_sf > 3.5);
if(~isempty(i))
  nthslice = 1;
  for slice = s.firstslice:s.lastslice
    nout_over(nthslice) = length(find(i==(slice+1)));
    nthslice = nthslice + 1;
  end
end
twf_over_alpha = (inv(X'*X)*X'*twf_over')';
twf_over_trend = twf_over_alpha(:,2);

% -------- Compute global stats -------------- %
glb_nover = sum(nover);
glb_twf_over = sum(twf_over .* repmat(nover, [1 s.ntp]))/glb_nover;
glb_twf_over_min = min(glb_twf_over);
glb_twf_over_max = max(glb_twf_over);
glb_twf_over_min_sf = sum(twf_over_min .* nover)/glb_nover;
glb_twf_over_max_sf = sum(twf_over_max .* nover)/glb_nover;
glb_twf_over_mean   = sum(twf_over_mean .* nover)/glb_nover;
glb_twf_over_stddev = std(glb_twf_over);
glb_twf_over_absdev = abs(glb_twf_over - glb_twf_over_mean);
glb_twf_over_z = glb_twf_over_absdev/glb_twf_over_stddev;
[m i] = max(glb_twf_over_z);
glb_twf_over_z_max  = m;
glb_twf_over_iz_max = i;
glb_twf_over_meanabsdev = sum(twf_over_meanabsdev .* nover)/glb_nover;
glb_twf_over_stddev_sf = glb_twf_over_meanabsdev;
glb_twf_over_snr_sf = glb_twf_over_mean./glb_twf_over_meanabsdev;
glb_nout_over = sum(nout_over);
glb_twf_over_alpha = (inv(X'*X)*X'*glb_twf_over')';
glb_twf_over_trend = glb_twf_over_alpha(:,2);

%-------------  SubThreshold: -------------------------%
% --------Compute slice-based stats ------------- %
twf_under_min  = min(twf_under')';
twf_under_max  = max(twf_under')';
twf_under_mean = mean(twf_under')';
twf_under_stddev = std(twf_under')';
twf_under_absdev = abs(twf_under - repmat(twf_under_mean,[1 s.ntp]));
twf_under_meanabsdev = mean(twf_under_absdev')';
twf_under_stddev_sf = twf_under_meanabsdev;
tmp = twf_under_stddev;
ind0 = find(tmp == 0);
tmp(ind0) = 10^10;
twf_over_z      = twf_over_absdev./(repmat(tmp,[1 s.ntp]));
twf_under_z      = twf_under_absdev./(repmat(twf_under_stddev,[1 s.ntp]));
twf_under_z_max  = max(twf_under_z')';
twf_under_snr_sf = twf_under_mean./twf_under_meanabsdev;
twf_under_z_sf   = twf_under_absdev./(repmat(twf_under_meanabsdev,[1 s.ntp]));
nout_under = zeros(s.nslices,1);
[i j] = find(twf_under_z_sf > 3.5);
if(~isempty(i))
  nthslice = 1;
  for slice = s.firstslice:s.lastslice
    nout_under(nthslice) = length(find(i==(slice+1)));
    nthslice = nthslice + 1;
  end
end
twf_under_alpha = (inv(X'*X)*X'*twf_under')';
twf_under_trend = twf_under_alpha(:,2);

% -------- Compute global stats -------------- %
glb_nunder = sum(nunder);
glb_twf_under = sum(twf_under .* repmat(nunder, [1 s.ntp]))/glb_nunder;
glb_twf_under_min = min(glb_twf_under);
glb_twf_under_max = max(glb_twf_under);
glb_twf_under_min_sf = sum(twf_under_min .* nunder)/glb_nunder;
glb_twf_under_max_sf = sum(twf_under_max .* nunder)/glb_nunder;
glb_twf_under_mean = sum(twf_under_mean .* nunder)/glb_nunder;
glb_twf_under_stddev  = std(glb_twf_under);
glb_twf_under_absdev = abs(glb_twf_under - glb_twf_under_mean);
glb_twf_under_z = glb_twf_under_absdev/glb_twf_under_stddev;
[m i] = max(glb_twf_under_z);
glb_twf_under_z_max  = m;
glb_twf_under_iz_max = i;
glb_twf_under_meanabsdev = sum(twf_under_meanabsdev .* nunder)/glb_nunder;
glb_twf_under_stddev_sf = glb_twf_under_meanabsdev; 
glb_twf_under_snr_sf = glb_twf_under_mean./glb_twf_under_meanabsdev;
glb_nout_under = sum(nout_under);
glb_twf_under_alpha = (inv(X'*X)*X'*glb_twf_under')';
glb_twf_under_trend = glb_twf_under_alpha(:,2);

ntot = glb_nover + glb_nunder;
glb_mean = (glb_twf_over_mean * glb_nover + ...
            glb_twf_under_mean * glb_nunder)/ntot;


% ------------- Over and Under ------------------%
ROverUnder = glb_twf_over_mean/glb_twf_under_mean;

% ------------ Over/Under Correlation ------- %
% Convert twf into a Unit-Normal Vector
a = glb_twf_over;
an0 = (a-mean(a));
glb_twf_over_norm = an0/sqrt(sum(an0.^2));
a = glb_twf_under;
an0 = (a-mean(a));
glb_twf_under_norm = an0/sqrt(sum(an0.^2));
% Correlation %
CorOverUnder = glb_twf_over_norm*glb_twf_under_norm'; %'
% Signficance %
a = glb_twf_over_norm;
b = glb_twf_under_norm;
eCor = (a'*a-eye(size(a'*a)))*b'; %'
eCorStd = std(eCor);
tCor = CorOverUnder/eCorStd;
tSigCor = tTest(s.ntp-2,tCor);


%------------- print meanval file -----------------%
fid = fopen(s.meanvalfile,'w');
if(fid == -1)
  fprintf('ERROR: could not open %s for writing\n',s.meanvalfile);
  return;
end
fprintf(fid,'%f\n',glb_twf_over_mean);
fclose(fid);

nvoxs = s.nrows * s.ncols * s.nslices;

%------------- print report -------------------- %
% -------- print report --------------- %
fid = fopen(s.reportfile,'w');
fprintf(fid,'# FS-FAST Intensity Normalization Report\n');
fprintf(fid,'# version: %s\n',version);
fprintf(fid,'# date:    %s\n',date);
fprintf(fid,'# Input Volume       %s\n',s.invol);
fprintf(fid,'# nrows              %d\n',s.nrows);
fprintf(fid,'# ncols              %d\n',s.ncols);
fprintf(fid,'# nslices            %d\n',s.nslices);
fprintf(fid,'# ntrs               %d\n',s.ntp);
fprintf(fid,'# inplaneres         %f\n',s.inplaneres);
fprintf(fid,'# betplaneres        %f\n',s.betplaneres);
fprintf(fid,'# TR                 %f\n',s.TR);
fprintf(fid,'# seqname            %s\n',s.seqname);
fprintf(fid,'# \n');

fprintf(fid,'# GlobalMean         %f\n',globalmean);
fprintf(fid,'# Relative Threshold Over %f\n',s.relthresh);
fprintf(fid,'# Absolute Threshold Over %f\n',s.absthresh);
fprintf(fid,'# Relative Threshold Under %f\n',s.relthresh_under);
fprintf(fid,'# Absolute Threshold Under %f\n',s.absthresh_under);
fprintf(fid,'# \n');

fprintf(fid,'# Over-Threshold Stats\n');
fprintf(fid,'# OV NVox        %d\n', glb_nover);
fprintf(fid,'# OV PctVox      %6.2f\n', 100*glb_nover/nvoxs);
fprintf(fid,'# OV Mean        %f\n', glb_twf_over_mean);
fprintf(fid,'# OV StdDev      %f\n', glb_twf_over_stddev);
fprintf(fid,'# OV AvgAbsDev   %f\n', glb_twf_over_stddev_sf);
fprintf(fid,'# OV Min         %f\n', glb_twf_over_min);
fprintf(fid,'# OV Max         %f\n', glb_twf_over_max);
fprintf(fid,'# OV Range       %f\n', glb_twf_over_max-glb_twf_over_min);
fprintf(fid,'# OV SNR         %f\n', glb_twf_over_mean/glb_twf_over_stddev);
fprintf(fid,'# OV ZAvg        %f\n',glb_twf_over_stddev_sf/glb_twf_over_stddev);
fprintf(fid,'# OV ZMax        %f\n', glb_twf_over_z_max);
fprintf(fid,'# OV ZMax Index  %d\n', glb_twf_over_iz_max);
fprintf(fid,'# OV Drift       %f\n', glb_twf_over_trend);
fprintf(fid,'# \n');

fprintf(fid,'# Under-Threshold Stats\n');
fprintf(fid,'# UN NVox        %d\n', glb_nunder);
fprintf(fid,'# UN PctVox      %6.2f\n', 100*glb_nunder/nvoxs);
fprintf(fid,'# UN Mean        %f\n', glb_twf_under_mean);
fprintf(fid,'# UN StdDev      %f\n', glb_twf_under_stddev);
fprintf(fid,'# UN AvgAbsDev   %f\n', glb_twf_under_stddev_sf);
fprintf(fid,'# UN Min         %f\n', glb_twf_under_min);
fprintf(fid,'# UN Max         %f\n', glb_twf_under_max);
fprintf(fid,'# UN Range       %f\n', glb_twf_under_max-glb_twf_under_min);
fprintf(fid,'# UN SNR         %f\n', glb_twf_under_mean/glb_twf_under_stddev);
fprintf(fid,'# UN ZAvg        %f\n', glb_twf_under_stddev_sf/glb_twf_under_stddev);
fprintf(fid,'# UN ZMax        %f\n', glb_twf_under_z_max);
fprintf(fid,'# UN ZMax Index  %d\n', glb_twf_under_iz_max);
fprintf(fid,'# UN Drift       %f\n', glb_twf_under_trend);
fprintf(fid,'# \n');

fprintf(fid,'# Over/Under Stats f\n');
fprintf(fid,'# OU Mean          %f\n', ROverUnder);
fprintf(fid,'# OU Cor           %f\n', CorOverUnder);
fprintf(fid,'# OU eCorStd       %f\n', eCorStd);
fprintf(fid,'# OU tCor          %f\n', tCor);
fprintf(fid,'# OU tSigCor       %f\n', tSigCor);
fprintf(fid,'# OU log10tSigCor  %f\n', -log10(tSigCor));
fprintf(fid,'# \n');

fprintf(fid,'# PctUnaccounted %6.2f\n',100*(1-(glb_nunder+glb_nover)/nvoxs));
fprintf(fid,'# \n');

fprintf(fid,'# StackFix-Based Stats\n');
fprintf(fid,'# SF Mean      %f\n', glb_twf_over_mean);
fprintf(fid,'# SF StdDev    %f\n', glb_twf_over_stddev_sf);
fprintf(fid,'# SF SNR       %f\n', glb_twf_over_snr_sf);
fprintf(fid,'# SF Min       %f\n', glb_twf_over_min_sf);
fprintf(fid,'# SF Max       %f\n', glb_twf_over_max_sf);
fprintf(fid,'# SF NOut      %d/%d\n', glb_nout_over,s.nslices*s.ntp);
fprintf(fid,'# \n');
fprintf(fid,'##Slc NOver  Mean    SFStd  SFSNR    Min     Max SFOut  ZMax  Trend\n');

slices = [s.firstslice:s.lastslice]'; %'
tmp = [slices nover twf_over_mean twf_over_stddev_sf twf_over_snr_sf ...
       twf_over_min twf_over_max nout_over twf_over_z_max twf_over_trend];

fprintf(fid,'%3d   %4d  %6.2f %7.4f %6.2f  %6.2f  %6.2f %2d  %6.2f %7.4f\n',tmp');%'

fclose(fid);

% Print some stuff to stdout %
fprintf(1,'# Global Mean        %f\n',globalmean);
fprintf(1,'# Mean Over/Under    %f\n', ROverUnder);
fprintf(1,'# Cor  Over/Under    %f\n', CorOverUnder);
fprintf(1,'# tSigCor Over/Under %f (%f)\n', tSigCor,-log10(tSigCor));

if(~isempty(s.twfstem))

  t = [0:s.ntp-1];

  fname = sprintf('%s-over',s.twfstem);
  fid = fopen(fname,'w');
  tmp = [t; glb_twf_over_norm; glb_twf_over; twf_over];
  fmt = repmat('%6.2f ',[1 size(tmp,1)]);
  fmt = [fmt '\n'];
  fprintf(fid,fmt,tmp);
  fclose(fid);

  fname = sprintf('%s-under',s.twfstem);
  fid = fopen(fname,'w');
  tmp = [t; glb_twf_under_norm; twf_under];
  fmt = repmat('%6.2f ',[1 size(tmp,1)]);
  fmt = [fmt '\n'];
  fprintf(fid,fmt,tmp);
  fclose(fid);

  % fname = sprintf('%s-overunder',s.twfstem);
  % fid = fopen(fname,'w');
  % tmp1 = glb_twf_over_norm;
  % tmp2 = glb_twf_under_norm;
  % tmp = [t; tmp1; tmp2];
  % fmt = repmat('%6.2f ',[1 size(tmp,1)]);
  % fmt = [fmt '\n'];
  % fprintf(fid,fmt,tmp);
  % fclose(fid);

end

r = 0;
fprintf(1,'fast_inorm: completed SUCCESSFULLY (%g)\n',toc);

return;
%----------------------------------------------------------%
%----------------------------------------------------------%
%----------------------------------------------------------%
%----------------------------------------------------------%


%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
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

      case '-i',
        arg1check(flag,narg,ninputargs);
        s.invol = inputargs{narg};
        narg = narg + 1;

      case {'-report','-reportfile'},
        arg1check(flag,narg,ninputargs);
        s.reportfile = inputargs{narg};
        narg = narg + 1;

      case {'-twf'},
        arg1check(flag,narg,ninputargs);
        s.twfstem = inputargs{narg};
        narg = narg + 1;

      case {'-meanval','-meanvalfile'},
        arg1check(flag,narg,ninputargs);
        s.meanvalfile = inputargs{narg};
        narg = narg + 1;

      case {'-firstslice', '-fs'}
        arg1check(flag,narg,ninputargs);
        s.firstslice = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-nslices', '-ns'}
        arg1check(flag,narg,ninputargs);
        s.nslices = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-TR'}
        arg1check(flag,narg,ninputargs);
        s.TR = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-inplaneres','-ps','-ipr'}
        arg1check(flag,narg,ninputargs);
        s.inplaneres = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-betplaneres','-st','-bpr'}
        arg1check(flag,narg,ninputargs);
        s.betplaneres = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-seqname'}
        arg1check(flag,narg,ninputargs);
        s.seqname = inputargs{narg};
        narg = narg + 1;

      case {'-relthresh','-thresh'}
        arg1check(flag,narg,ninputargs);
        s.relthresh = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-absthresh'}
        arg1check(flag,narg,ninputargs);
        s.absthresh = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case '-o',
        arg1check(flag,narg,ninputargs);
        s.outvol = inputargs{narg};
        narg = narg + 1;

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
  fprintf(1,'  fast_inorm\n');
  fprintf(1,'     -i  invol       : volume to normalize\n');
  fprintf(1,'     -TR val         : TR \n');
  fprintf(1,'     -inplaneres  val : pixel size\n');
  fprintf(1,'     -betplaneres val : between plane resolution\n');
  fprintf(1,'     -seqname  name  : name of acquisition sequence\n');
  fprintf(1,'     -relthresh val  : relative threshold (default .75) \n');
  fprintf(1,'     -absthresh val  : absolute threshold \n');
  fprintf(1,'     -report filename   : name of file for report\n');
  fprintf(1,'     -meannval filename : name of file in which to store mean value\n');
  %fprintf(1,'     -o  outvol     : \n');
  %fprintf(1,'     -rescale   val : target value (defaut 1000)\n');

return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = inorm_struct
  s.invol       = '';
  s.involendian = '';
  s.involext    = '';
  s.firstslice  = 0;
  s.lastslice   = 0;
  s.nslices    = -1;
  s.nrows      = 0;
  s.ncols      = 0;
  s.ntp      = 0;
  s.outvol     = '';
  s.relthresh  = [];
  s.absthresh  = [];
  s.relthresh_under  = 0.25;
  s.absthresh_under  = [];
  s.rescale    = 1000;
  s.reportfile  = [];
  s.meanvalfile = [];
  s.twfstem     = [];
  s.verbose     = 0;
  s.TR          = -1;
  s.inplaneres   = -1;
  s.betplaneres  = -1;
  s.seqname     = 'unknown';
return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Print data structure
function s = print_struct(s,fid)
  if(nargin == 1) fid = 1; end

  fprintf(fid,'invol        %s\n',s.invol);
  fprintf(fid,'firstslice   %d\n',s.firstslice);
  fprintf(fid,'nslices      %d\n',s.nslices);
  %fprintf(fid,'outvol       %s\n',s.invol);
  fprintf(fid,'relthresh    %f\n',s.relthresh);
  if(~isempty(s.absthresh))
    fprintf(fid,'absthresh    %f\n',s.absthresh);
  end
  %if(~isempty(s.rescale))
  %  fprintf(fid,'rescale      %f\n',s.rescale);
  %end
  fprintf(fid,'reportfile   %s\n',s.reportfile);
  fprintf(fid,'meanvalfile  %s\n',s.meanvalfile);
  fprintf(fid,'twfstem      %s\n',s.twfstem);
  fprintf(fid,'TR           %f\n',s.TR);

return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Check argument consistency, etc %%%
function s = check_params(s)

  fprintf(1,'Checking Parameters\n');
  if(size(s.invol,1) < 1) 
    msg = sprintf('ERROR: must have an input volume');
    qoe(msg);error(msg);
  end
  [s.nrows s.ncols s.ntp fs ns s.involendian s.involext] = ...
    fmri_bfiledim(s.invol);
  if(s.nslices == -1)      s.nslices = ns; end
  s.lastslice = s.firstslice + s.nslices - 1;
  %s.lastplane = s.firstplane + s.nplanes - 1;

  if(s.nslices < 1) 
     msg = sprintf('ERROR: nslices = %d, must be > 0',s.nslices);
     qoe(msg);error(msg);
  end
  if(s.firstslice < 0) 
     msg = sprintf('ERROR: firstslice (%d) < 0',s.firstslice);
     qoe(msg);error(msg);
  end
  if(s.lastslice >= ns) 
     msg = sprintf('ERROR: last slices (%d) >= nslices avail (%d)',...
                   s.lastslice,ns);
     qoe(msg);error(msg);
  end

  if(s.TR <= 0) 
    msg = sprintf('ERROR: TR = %f, must be > 0',s.TR);
    qoe(msg);error(msg);
  end

  if(strcmp(s.invol,s.outvol))
    if(s.verbose) fprintf(1,'INFO: overwriting input volume\n'); end
  end

  if(isempty(s.relthresh) & isempty(s.absthresh))
     s.relthresh = .75;
  end

  if(isempty(s.reportfile))
     s.reportfile = sprintf('%s.report',s.invol);
  end

  if(isempty(s.meanvalfile))
     s.meanvalfile = sprintf('%s.meanval',s.invol);
  end

  if(isempty(s.twfstem))
    s.twfstem = sprintf('%s.twf',s.invol);
  end

return;


