% HDRFiles, datFiles, OutFile
% Effect, Weight


%
% fmri_axaslice.m
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

if( ~exist('QuitOnError') ) QuitOnError = 1; end

%%% ---- Check that all the variables are defined --- %%%
VarNameList = [];
VarNameList = strvcat(VarNameList,'HDRFiles');
VarNameList = strvcat(VarNameList,'OutputFile');
VarNameList = strvcat(VarNameList,'OutputDatFile');
VarNameList = strvcat(VarNameList,'datFile');
VarNameList = strvcat(VarNameList,'Effect');
VarNameList = strvcat(VarNameList,'Weight');
nVar = size(VarNameList,1);
for n = 1:nVar,
  if( exist(deblank(VarNameList(n,:))) ~= 1)
    msg = sprintf('Error: Variable %s does not exist',VarNameList(n,:));
    qoe(msg);error(msg);
  end
end


Effect = lower(Effect);
Weight = lower(Weight);

if(~strcmp(Effect,'random') & ~strcmp(Effect,'fixed'))
  msg = sprintf('Effect %s unkown, should be random or fixed');
  qoe(msg);error(msg);
end

if(~strcmp(Weight,'subject') & ~strcmp(Weight,'event'))
  msg = sprintf('Weight %s unkown, should be subject or event');
  qoe(msg);error(msg);
end


%%% -------- read in the dat file ---------- %%
[nNNCond nHEst DOF  TR nRuns nTP nRows nCols ...
   nSkip DTOrder Rescale TW tPreStim HanRad BASeg ...
   GammaFit gfDelta gfTau NullCondId SumXtX] =   fmri_lddat(datFile(1,:));

nCond = nNNCond + 1;
Nch = nNNCond*nHEst;
Ns = size(HDRFiles,1);

fprintf(1,'\nNs = %d, nCond = %d, Nch = %d, nHEst = %d\n',...
        Ns, nCond, Nch, nHEst);

hselavg = fmri_ldbfile(HDRFiles);
nRows  = size(hselavg,1);
nCols  = size(hselavg,2);

havg = zeros(nRows,nCols,Nch,Ns);
hstd = zeros(nRows,nCols,Nch,Ns);

eVar = 0;
for s = 1:Ns,
  fprintf(1,' --------------- s = %d ------------\n',s);
  [ha eVar2 hs] = fmri_sa2sxa(hselavg(:,:,:,s),nHEst);
  havg(:,:,:,s) = ha;
  hstd(:,:,:,s) = hs;
  eVar2 = eVar2 + eVar;
  fprintf(1,'havg: %g %g\n', min(reshape1d(ha)),max(reshape1d(ha)));
  fprintf(1,'hstd: %g %g\n', min(reshape1d(hs)),max(reshape1d(hs)));
end

if(strcmp(Weight,'subject'))
  [haa hsa DOF Ms] = fmri_avgxavg(Effect,havg,hstd);
else
  [haa hsa DOF Ms] = fmri_avgxavg(Effect,havg,hstd,SumXtX);
end

fprintf(1,' --------------- pooled ------------\n');
fprintf(1,'havg: %g %g\n',min(reshape1d(haa)),max(reshape1d(haa)));
fprintf(1,'hstd: %g %g\n',min(reshape1d(hsa)),max(reshape1d(hsa)));

Nch2 = nCond * nHEst;

haa = reshape(haa, [nRows nCols nNNCond nHEst]);
hsa = reshape(hsa, [nRows nCols nNNCond nHEst]);

haa = permute(haa, [4 3 2 1]);
hsa = permute(hsa, [4 3 2 1]);

ysa = zeros(nRows, nCols, nCond, 2, nHEst);
ysa = permute(ysa, [5 3 2 1 4]);
ysa(:,[2:nCond],:,:,1) = haa;
ysa(:,[2:nCond],:,:,2) = hsa;
ysa = permute(ysa, [1 5 2 3 4]);
ysa = reshape(ysa, [2*Nch2 nRows nCols]);
ysa = permute(ysa, [2 3 1]);

fmri_svbfile(ysa,OutputFile);
nNoiseAC = -1;
fmri_svdat(OutputDatFile,TR,TW,tPreStim,nNNCond,DOF,...
           nRuns,nTP,nRows,nCols,nSkip,DTOrder,Rescale,HanRad,...
           nNoiseAC,BASeg,GammaFit,gfDelta,gfTau,NullCondId,SumXtX);
return;


fid = fopen(OutputDatFile,'w');
if(fid == -1)
  msg = sprintf('Could not open %s for writing',OutputDatFile);
  qoe(msg);error(msg);
end
fprintf(fid,'avgxavg dat file\n');
fprintf(fid,'TR         %g\n',TR);
fprintf(fid,'TimeWindow %g\n',TW);
fprintf(fid,'TPreStim   %g\n',tPreStim);
fprintf(fid,'nCond      %d\n',nNNCond+1);
fprintf(fid,'nHEst      %d\n',nHEst);
fprintf(fid,'DOF        %d\n',DOF);
fprintf(fid,'nSessions  %d\n',Ns);
fprintf(fid,'nRows      %d\n',nRows);
fprintf(fid,'nCols      %d\n',nCols);
fclose(fid);

return;

r = 43;
c = 42;
ha1 = squeeze(haa(r,c,1,1:nHEst));
hs1 = squeeze(hsa(r,c,1,1:nHEst));
ha2 = squeeze(haa(r,c,1,nHEst+1:2*nHEst));
hs2 = squeeze(hsa(r,c,1,nHEst+1:2*nHEst));
figure(1);
errorbar(TR*[0:nHEst-1],ha1,hs1);
hold
errorbar(TR*[0:nHEst-1],ha2,hs2);
hold


