% fmri_roiavg
%
% 

%
% fmri_roiavg.m
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


%---------------------------------------------- %
if(0)
  TopDir  =   '/space/raid/5/users/caplan/sujith'
  hdrstem = strcat(TopDir,'/990215MG/selavg-dt/talairach2/th');
  sigstem = strcat(TopDir,'/average/allh15dt-fe/stats/pm-allfix-mdiff');
  threshold = .000001;
  sigformat = 'ln';
  roi_rowcol = [12 6 18 12];
  FirstSlice =  19;
  nSlices    =   3;
  sigframe   =  1;
  sigsign    = '+';

  %TopDir  =   '/space/raid/5/users/greve/synthtst/AllSig/selavg/'
  %hdrstem = strcat(TopDir,'h01');
  %roi_rowcol = [12 6 18 12];
  %FirstSlice =  0;
  %nSlices    =  1;
  ShowResults = 1;
  OutStem = 'roi1avg';
end
%---------------------------------------------- %

fprintf('\n\n%s\n',hdrstem);

havgall = [];
hstdall = [];
NvPerSlice = [];
switch(sigformat)
  case 'ln',    sigthreshold = -log(threshold);
  case 'log10', sigthreshold = -log10(threshold);
  case 'raw',   sigthreshold = 1-threshold;
  otherwise,
    msg = sprintf('Sig Format: %s not recognized (use: ln, log10, or raw)');
    qoe(msg);error(msg);
end

fprintf('Sig Threshold for %s is %g\n',sigformat,sigthreshold);

Nvmax = 0;
for slc = FirstSlice : FirstSlice+nSlices-1,
  fprintf(1,'Slice  %2d  --------------- \n',slc);

  datfile = sprintf('%s.dat',hdrstem);
  hdrdat = fmri_lddat3(datfile);

  hsafile = sprintf('%s_%03d.bfloat',hdrstem,slc);
  hsa = fmri_ldbfile(hsafile);
  [havg hstd hdrdat] = fmri_untangle(hsa,hdrdat);
  %havg = randn(hdrdat.Nrows,hdrdat.Ncols,Nch);
  %havglist(:,:,:,s) = havg;

  Nv = hdrdat.Nrows * hdrdat.Ncols;
  Nch = hdrdat.Nc * hdrdat.Nh;
  havg = permute(havg, [4 3 1 2]);
  havg = reshape(havg, [Nch Nv]);

  hstd = permute(hstd, [4 3 1 2]);
  hstd = reshape(hstd, [Nch Nv]);

  roi_index = roi2ind([hdrdat.Nrows hdrdat.Ncols],roi_rowcol);
  fprintf('  Number of indicies in ROI %3d %3d %3d %3d:  %d\n',...
          roi_rowcol,length(roi_index));
  Nvmax = Nvmax + length(roi_index);

  if(~isempty(sigstem) & abs(threshold) < 1)
    roimask = zeros(hdrdat.Nrows,hdrdat.Ncols);
    roimask(roi_index) = 1;
    sigfile = sprintf('%s_%03d.bfloat',sigstem,slc);
    %fprintf('   Loading %s\n',sigfile);
    p = fmri_ldbfile(sigfile);
    if(sigframe > size(p,3))
      msg = sprintf('sigframe = %d > sigvol depth = %d\n',sigframe, size(p,3));
      qoe(msg);error(msg);
    end
    p = squeeze(p(:,:,sigframe)).*roimask;
    switch(sigsign)
      case '+',   throi_index = find(p      >  sigthreshold);
      case '-',   throi_index = find(p      < -sigthreshold);
      case '+/-', throi_index = find(abs(p) >  sigthreshold);
      otherwise
        msg = sprintf('Sig Sign: %s not recognized (use: +, -, or +/-)');
      qoe(msg);error(msg);
    end
    fprintf('  Number of voxels in ROI above threshold  %d\n',...
          length(throi_index));
    roi_index = throi_index;
  end

  havgall = [havgall havg(:,roi_index)];
  hstdall = [hstdall hstd(:,roi_index)];
  NvPerSlice = [NvPerSlice length(roi_index)];

end % for slc

dof1 = size(havgall,2);
fprintf('Nv over threshold: %d/%d\n',dof1,Nvmax);
dof  = dof1*ones(hdrdat.Nc,1);
if(dof1==0)
  hAvg = zeros([hdrdat.Nc hdrdat.Nh]);	
  hStd = zeros([hdrdat.Nc hdrdat.Nh]);	
else
  hAvg = mean(havgall,2);
  hAvg = reshape(hAvg, [hdrdat.Nh hdrdat.Nc])';	%'
  hStd = sqrt(mean(hstdall.^2,2));
  %hStd = std(havgall,[],2);
  hStd = reshape(hStd, [hdrdat.Nh hdrdat.Nc])';	%'
  if(white) hStd = hStd/sqrt(dof1); end
end

t = hdrdat.TR*[0:hdrdat.Nh-1] - hdrdat.TPreStim;

if(ShowResults)
  fprintf('showing results\n');
  hHDR = figure(1);
  hdrviewlite('init',hHDR,t,dof);
  hdrviewlite('plot',hHDR,hAvg,hStd, [0 0]);
  uiwait(hHDR);
  fprintf('done showing results\n');
end

fprintf('OutStem %s\n',OutStem);
save(OutStem,'hdrstem','t','dof','hAvg','hStd',...
     'hdrdat','roi_rowcol','FirstSlice','nSlices','threshold',...
     'sigthreshold','sigsign','NvPerSlice','sigformat','Nv','Nvmax');

if(report)
  h0 = repmat(hAvg(1,:), [hdrdat.Nc 1]);
  htmp = hAvg-h0;
  htmp = htmp(2:hdrdat.Nc,:);
  RepFile = sprintf('%s.txt',OutStem);
  fmt = repmat('%8.4f ',[1 hdrdat.Nh]);
  fmt = [fmt '\n'];
  fid = fopen(RepFile,'w');
  fprintf(fid,fmt,htmp'); %'
  fclose(fid);
end


fprintf('fmri_roiavg: done\n\n');
if(QuitOnError) quit; end;

