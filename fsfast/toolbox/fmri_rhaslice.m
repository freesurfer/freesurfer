% fmri_rhaslice.m
% Computes the average HDR of a ROI
%
% Variables:
% tAvgFile, IndexFile, HDRFile, HDRDatFile, OutputFile, Init
%
%


%
% fmri_rhaslice.m
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

fprintf('-------- fmri_rhaslice.m ------------\n');
fprintf('HDRFile = %s\n',HDRFile);
fprintf('Init = %d\n',Init);

%% Load the HDR File %%%
hdr = fmri_ldbfile(HDRFile);
[nRows nCols Nch] = size(hdr);
nVoxels = nRows*nCols;
hdr = reshape(hdr,[nVoxels Nch])';

%% Load the tAvg File %%%
tavg = fmri_ldbfile(tAvgFile);
tavg = reshape1d(tavg);

%% Load the datfile for the HDR File
fid=fopen(deblank(HDRDatFile),'r');
if( fid == -1 )
  msg = sprintf('Could not open dat file %s\n',HDRDatFile);
  qoe(msg);  error(msg);
end
fscanf(fid,'%s',1);
TR      = fscanf(fid,'%f',1);
fscanf(fid,'%s',1);
TW      = fscanf(fid,'%f',1);
fscanf(fid,'%s',1);
tPreStim     = fscanf(fid,'%f',1);
fscanf(fid,'%s',1);
nCond   = fscanf(fid,'%d',1);
fscanf(fid,'%s',1);
nHEst   = fscanf(fid,'%d',1); % nperevent
fclose(fid);
nPreStim     = floor(tPreStim/TR);
nNNC    = nCond - 1;

%% Load the index file %%
fid = fopen(deblank(IndexFile),'r');
if(fid == -1)
  msg = sprintf('Could not open %s',IndexFile);
  qoe(msg);error(msg);
end
[ind nind] = fscanf(fid,'%d');
fclose(fid);
if(mod(nind,2) ~= 0)
  msg = sprintf('Odd number of elements in index file (%d)',nind);
  qoe(msg);error(msg);
end
ind = ceil(ind/2); % convert from structural to functional
ind = reshape(ind, [2 nind/2])'; %'
ind = sub2ind([nRows nCols], ind(:,1), ind(:,2)); % convert from subscripts
nVoxComb = nind/2;

%% get averages of temporal averages %%
tavgavg = mean(tavg(ind));

%% get the hdrs for the relevant voxels %%
hdr = hdr(:,ind);

%% reshape into convenient dimensions %%
hdr = reshape(hdr,[nHEst 2 nCond nVoxComb]);  

%% compute difference with condition 0 %%
h = [];
for c = 1:nNNC,
  h(:,c,:) = hdr(:,1,c+1,:) - hdr(:,1,1,:);
end

%% compute averages and std across voxels %%
havg = mean(h,3);
% hstd = std(h,[],3);

%% Check whether to init or accumulate %%
if(~Init) 
  fid = fopen(OutputFile,'r');
  if(fid==-1)
    msg = sprintf('Cannot open %s for reading',OutputFile);
    qoe(msg);error(msg);
  end

  q = fscanf(fid,'%f');
  fclose(fid);
  q = reshape(q, [3+nNNC nHEst])';
  nVoxCombOld = q(1,2);
  tavgavgOld  = q(1,3);
  havgOld = q(:,4:4+nNNC-1);

  hsum = havg*nVoxComb + havgOld*nVoxCombOld;
  tsum = tavgavg*nVoxComb + tavgavgOld*nVoxCombOld;

  nVoxComb = nVoxComb + nVoxCombOld;
  havg = hsum/nVoxComb;
  tavgavg = tsum/nVoxComb;
end


%%% Save results %%%%%%%%%%%%%%%
fid = fopen(OutputFile,'w');
if(fid==-1)
  msg = sprintf('Cannot open %s for writing',OutputFile);
  qoe(msg);error(msg);
end

%% Format: t, DOF, tavgavg, avgsig_c1 avgsig_c2 ... %%
for n = 1:nHEst
  t = TR*(n-1);
  fprintf(fid,'%5.1f %6d %6.1f ',t,nVoxComb,tavgavg);
  fprintf(fid,'%5.2f ',havg(n,:));

%  fprintf(fid,'%8.5f ',havg(n,:)/tavgavg);
%  for c = 1:nNNC,
%    fprintf(fid,'%5.2f ',hstd(n,c));
%  end

  fprintf(fid,'\n');
end

fclose(fid);


