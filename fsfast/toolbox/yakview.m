% yakview - views images, stat overlays, and hemodynamic responses.
% $Id: yakview.m,v 1.5 2004/03/22 20:29:30 greve Exp $

fprintf(1,'\n\n');
fprintf(1,'yakview: $Id: yakview.m,v 1.5 2004/03/22 20:29:30 greve Exp $\n');

if(~exist('UseVersion')) UseVersion = 2; end

yakTitle = ImgFile;

ln10 = log(10.0); % natural log of 10.0
p = [];
Ch = [];

% Make sure the min threshold is not greater than the max
if(SigThresh > SigThreshMax) 
  msg = sprintf('SigThresh (%g) > SigThreshMax (%g) ',SigThresh,SigThreshMax);
  qoe(msg); error(mag);
end

fprintf(1,'\n\n');
fprintf(1,'Base %s \n',ImgFile);
fprintf(1,'Overlay %s \n',SigFile);
fprintf(1,'Hemodynamic Data %s \n',HDRFile);
fprintf(1,'Offset %s \n',OffFile); 
fprintf(1,'Raw Data %s \n',RawFile); 
fprintf(1,'RawFunc: %s', rawfunc);
fprintf(1,'\n\n');

fprintf(1,'Loading Base ...       '); tic;
if(~ImgMkMosaic)
  s = fast_ldbfile(ImgFile,1);
  if(isempty(s))
    fprintf('ERROR: loading %s\n',ImgFile);
    return;
  end
  
  fprintf('%g sec\n',toc);
  s = s(:,:,1); % keep only first plane of the base 
  basesize = size(s);
else
  if(fmtimg) s = fast_ldanalyze(ImgFile);
  else       s = fmri_ldbvolume(ImgFile);
  end
  fprintf('%g  sec\n',toc);
  s = s(:,:,:,1);
  s = permute(s, [2 3 1]);   % row, col, slice
  fprintf(1,'Making Base Mosaic ... '); tic;
  switch(MosaicDirection)
     case 'row'
       s = permute(s, [3 2 1]);  
     case 'col'
       s = permute(s, [3 1 2]);  
  end
  basesize = size(s);
  s = vol2mos(s);
  fprintf('%g  sec\n',toc);
end

if(ImgHistEQ)
  HistEQThr = .99;
  fprintf('Equalizing Base Image with threshold %g ... ',HistEQThr); tic;
  [s xthresh nclip psqueeze] = drsqueeze(s,HistEQThr);
  fprintf('%g  sec\n',toc);
  %fprintf('  xthresh  = %g\n',xthresh);
  %fprintf('  nclip    = %d / %d\n',nclip,prod(size(ymos(:,:,n))));
  %fprintf('  psqueeze = %g\n',psqueeze);
end

if(~isempty(SigFile))
  fprintf(1,'Loading Overlay\n'); tic;
  if(~SigMkMosaic)
    p = fmri_ldbfile(SigFile);
    if(pneg) p = -p; end
    if(~isempty(SigMaskFile))
      fprintf(1,'Loading Mask \n'); tic;
      pmask = fmri_ldbfile(SigMaskFile);
      pmask = abs(pmask) > SigMaskThresh;
      pmask = repmat(pmask,[1 1 size(p,3)]);
      p = p.*pmask;
      clear pmask
    end
  else
    if(fmtimg) p = fast_ldanalyze(SigFile);
    else       p = fmri_ldbvolume(SigFile);
    end
    if(~isempty(cutends))  p([1 size(p,1)],:,:) = cutends;  end
    if(~isempty(SigMaskFile))
      if(fmtimg) pmask = fast_ldanalyze(SigMaskFile);
      else       pmask = fmri_ldbvolume(SigMaskFile);
      end
      pmask = abs(pmask) > SigMaskThresh;
      pmask = repmat(pmask,[1 1 1 size(p,4)]);
      p = p.*pmask;
      clear pmask
    end
    fprintf('%g  sec\n',toc);
    fprintf(1,'Making Overlay Mosaic ... '); tic;
    p = permute(p, [2 3 1 4]);  
    switch(MosaicDirection)
       case 'row'
         p = permute(p, [3 2 1 4]);  
       case 'col'
         p = permute(p, [3 1 2 4]);  
    end
    p = vol2mos(p);
    fprintf('%g  sec\n',toc);
  end

  if(strcmp(SigFormat,'ln'))
    fprintf('Converting from ln to log10\n');
    p = -sign(p) .* log10( exp(-abs(p)));
  end

  yakTitle = strcat('SigMap:  ',SigFile,'  Struct:  ',ImgFile);
end


if(~exist('TR'))TR = 0;end

if(~isempty(datFile))  hd = fmri_lddat3(datFile); end
hsa = [];
hoffset = [];
Nnnc = 0;
TPS = 0;
dof = 0;
hd = [];
if(~isempty(HDRFile))
  fprintf(1,'Loading Hemodynamic Data ...      '); tic;
  if(~HDRMkMosaic)
    hsa = fmri_ldbfile(HDRFile);
    fprintf('%g sec\n',toc);
  else
    hsa = fmri_ldbvolume(HDRFile);
    fprintf('%g  sec\n',toc);
    fprintf(1,'Making Hemodynamic Data Mosaic ... '); tic;
    hsa = permute(hsa, [2 3 1 4]);  
    switch(MosaicDirection)
       case 'row'
         hsa = permute(hsa, [3 2 1 4]);  
       case 'col'
         hsa = permute(hsa, [3 1 2 4]);  
    end
    hsa = vol2mos(hsa);
    fprintf('%g  sec\n',toc);
  end

  %fprintf(1,'Loading DatFile %s\n',datFile);
  hd = fmri_lddat3(datFile);
  Nnnc = hd.Nnnc;
  TR   = hd.TER;
  TPS  = hd.TPreStim;

  if(~isempty(hd.hCovMtx))
    Ch = hd.hCovMtx;
  elseif(~isempty(hd.SumXtX))
    Ch = hd.SumXtX;
  else
    Ch = [];
  end

  if(hd.Version == 0) 
    fprintf(1,'Loading DOFFile \n');
    dof = fmri_lddof(HDRFile);
    dof = dof(:,2);
  else
    d = diag(hd.SumXtX);
    nndof = d(1:hd.Nh:length(d))';%'
    dof = [hd.DOF-sum(nndof) nndof];
  end

  if(isempty(SigFile))
    yakTitle = strcat('HDR:',HDRFile,'Struct:',ImgFile);
  end

  if(~isempty(OffFile))
    fprintf(1,'Loading Offset  ...      '); tic;
    if(~OffMkMosaic)
      hoffset = fmri_ldbfile(OffFile);
      fprintf('%g sec\n',toc);
      hoffset = hoffset(:,:,1);
    else
      if(fmtimg) hoffset = fast_ldanalyze(OffFile);
      else       hoffset = fmri_ldbvolume(OffFile);
      end
      fprintf('%g  sec\n',toc);
      hoffset = hoffset(:,:,:,1);
      fprintf(1,'Making Offset Mosaic ... '); tic;
      hoffset = permute(hoffset, [2 3 1 ]);  
      switch(MosaicDirection)
         case 'row'
           hoffset = permute(hoffset, [3 2 1]);  
         case 'col'
           hoffset = permute(hoffset, [3 1 2]);  
      end
      hoffset = vol2mos(hoffset);
      fprintf('%g  sec\n',toc);
    end
  end
end

yraw = [];
if(~isempty(RawFile))
  fprintf(1,'Loading Raw Data ...       '); tic;
  if(~RawMkMosaic)
    yraw = fmri_ldbfile(RawFile);
    fprintf('%g sec\n',toc);
  else
    yraw = fmri_ldbvolume(RawFile);
    fprintf('%g  sec\n',toc);
    fprintf(1,'Making Raw Data Mosaic ... '); tic;
    yraw = permute(yraw, [2 3 1 4]);  
    switch(MosaicDirection)
       case 'row'
         yraw = permute(yraw, [3 2 1 4]);  
       case 'col'
         yraw = permute(yraw, [3 1 2 4]);  
    end
    yraw = vol2mos(yraw);
    fprintf('%g  sec\n',toc);
  end
end

yakTitle = strcat('Yakview -- ',yakTitle);

if(~isempty(yvtitle))
  yakTitle = yvtitle;
end

close all;

fprintf(1,'\n\n');

if(UseVersion == 1)
yak('init',s,'-Nnnc',Nnnc,'-h',hsa,'-dof',dof,...
    '-TR',TR,'-tPreStim',TPS,'-A',p,'-pmin',SigThresh,...
    '-pmax',SigThreshMax,'-raw',yraw,'-off',hoffset,'-Ch',Ch,...
    '-hdat', hd,'-nskip',Nskip,'-rawfunc',rawfunc);
end
if(UseVersion == 2)
yak2('init',s,'-Nnnc',Nnnc,'-h',hsa,'-dof',dof,...
    '-TR',TR,'-tPreStim',TPS,'-A',p,'-pmin',SigThresh,...
    '-pmax',SigThreshMax,'-raw',yraw,'-off',hoffset,'-Ch',Ch,...
    '-hdat', hd,'-nskip',Nskip,'-rawfunc',rawfunc,...
    '-slice',SliceNo,basesize,'-hdatfile',datFile);
end




set(gcf,'Name',yakTitle);
%title(yakTitle);

clear p hsa s;

return;

