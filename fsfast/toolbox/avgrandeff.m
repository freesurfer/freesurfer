% avgrandeff
%
% Average using a Random Effects Model
%
% InputStems, SubtCond0, FirstSlice, nSlices, OutStem
% RescaleAvg, ROI


%
% avgrandeff.m
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

LastSlice = FirstSlice + nSlices - 1;
nInputs = size(InputStems,1);

%fprintf('%d %d %d\n\n',FirstSlice,LastSlice,nSlices);

indROI = [];

for slice = FirstSlice:LastSlice,

  fprintf(1,'--------- Processing Slice %d -----------\n',slice);

  if(isempty(ROI) | slice == FirstSlice)
    hsum  = [];
    hsum2 = [];
    dof   = 0;
  end

  for n = 1:nInputs,

    InStem = deblank(InputStems(n,:));
    fprintf(1,'   Processing Input %s, Slice %d\n',InStem,slice);

    datFile = sprintf('%s.dat',InStem);
    % fprintf(1,'Loading DatFile %s\n',datFile);
    [Nnnc,Nh,DOF,TR,nRuns,nTP,nRows,nCols,nSkip,DTOrder,...
       Rescale,TW,TPS,HanRad,BASeg,GammaFit, gfDelta, gfTau, ...
       NullCondId, SumXtX] = fmri_lddat(datFile);
    fprintf('\n');

    Nc = Nnnc + 1;
    indhavg = [];
    for c=1:Nc,
      indhavg = [indhavg ([1:Nh] + 2*Nh*(c-1))];
    end

    hfile      = sprintf('%s_%03d.bfloat',InStem,slice);
    hsa        = fmri_ldbfile(hfile);
    [nRows nCols Nch] = size(hsa);
    Nch = Nch/2;
    Nv = nRows*nCols;

    havgslice  = hsa(:,:,indhavg);

    havgslice  = reshape(havgslice, [Nv Nch])';
    if(SubtCond0)
      h0 = havgslice(1:Nh,:);
      h0 = repmat(h0, [Nc 1]);
      havgslice = havgslice - h0;
    end

    if(0) 
      havgslice  = reshape(havgslice', [nRows nCols Nch]);
      HanFilter = fmri_hankernel(3.0);
      fprintf('Hanning Filter\n');
      havgslice = fmri_spatfilter(havgslice,HanFilter);
      havgslice  = reshape(havgslice, [Nv Nch])';
    end

    if(0) 
      fprintf('   Fixing Baseline\n');
      havgslice = reshape(havgslice, [Nh Nc Nv]);
      blc = mean(havgslice([1:3],:,:),1);
      blc = repmat(blc, [Nh 1 1]);
      havgslice = havgslice - blc;
    end

    if(~isempty(ROI) & slice == FirstSlice & n == 1 )
      indROI = roi2ind([nRows nCols],ROI);
      fprintf('Region of Interest Indicies\n');
      fmt = repmat('%4d ',[1 10]);
      fmt = strcat(fmt,'\n');
      fprintf(fmt,indROI);
      fprintf('\n\n');
    end

    [hsum hsum2 dof] = fmri_accrandeff(havgslice,hsum,hsum2,dof,...
                                     RescaleAvg,indROI);

  end % Loop over slices %

  %%%% --- Compute and Save HDR fImage --- %%%%%
  if(isempty(ROI))

    [havg hstd] = fmri_avgrandeff(hsum,hsum2,dof);

    %%% Save Average of Averages (and header) %%%%%%
    OutFile = sprintf('%s_%03d.bfloat',OutStem,slice);
    a = reshape(havg, [Nh Nc*Nv]);
    b = reshape(hstd, [Nh Nc*Nv]);
    c = [a; b];
    d = reshape(c, [2*Nch Nv]); 
    hsa = reshape(d', [nRows nCols 2*Nch]);
    fprintf(1,'  Saving to %s\n',OutFile);
    fmri_svbfile(hsa,OutFile);

    %%% Save DOF %%%
    dofFile = sprintf('%s_%03d.dof',OutStem,slice);
    vdof = dof * ones(Nc,1);
    mdof = [ [0:Nnnc]' vdof vdof-1];
    fid = fopen(dofFile,'w');
    if(fid == -1)
      msg = sprintf('Could not open %s for writing',dofFile);
      qoe(msg);error(msg);
    end
    fprintf(fid,'%d %d %d\n',mdof');
    fclose(fid);

    %%% Save Dat File %%%
    datFile = sprintf('%s.dat',OutStem);
    if(SumXtX ~= 0)
      fmri_svdat(datFile,TR,TW,TPS,Nnnc,dof,nRuns,nTP,...
          nRows,nCols,nSkip,DTOrder,Rescale,HanRad,0,BASeg,...
          GammaFit,gfDelta,gfTau,NullCondId,SumXtX);
    else
      fid = fopen(datFile,'w');
      if(fid == -1)
        msg = sprintf('Could not open %s for writing',datFile);
        qoe(msg);error(msg);
      end
      fprintf(fid,'TR  %g\n',TR);
      fprintf(fid,'TW  %g\n',TW);
      fprintf(fid,'TPS %g\n',TPS);
      fprintf(fid,'Nc  %d\n',Nc);
      fprintf(fid,'Nh  %g\n',Nh);
      fclose(fid);

    end

  end

end % Loop over input files %



%%% Save as individual curves %%%%%%%%
if(~isempty(ROI))

  [havg hstd] = fmri_avgrandeff(hsum,hsum2,dof);

  t = TR*[0:Nh-1]' -TPS;
  havg = reshape(havg, [Nh Nc]);
  hstd = reshape(hstd, [Nh Nc]);
  hstderr = hstd/sqrt(dof-1);

  m = [t havg hstd hstderr];
  fmt = repmat('%g ', [1 size(m,2)]);
  fmt = strcat(fmt,'\n');

  OutFile = sprintf('%shdr.roi',OutStem);
  fid = fopen(OutFile,'w');
  if(fid == -1)
      msg = sprintf('Could not open %s for writing',OutFile);
      qoe(msg);error(msg);
  end
  fprintf(fid,'# %d \n',dof);
  fprintf(fid,fmt,m');
  fclose(fid);

  hhdrview = figure(1);
  vdof = repmat(dof, [1 Nc]);
  hdrviewlite('init',gcf,t',vdof);  
  hdrviewlite('plot',gcf,havg',hstd',[1 1]);
  uiwait(hhdrview);
end

fprintf('\n\n');

return;
