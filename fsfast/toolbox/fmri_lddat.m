function [nNNC,nHEst,DOF,TR,nRuns,nTP,nRows,nCols,nSkip,DTOrder,...
          Rescale,TW,TPS,HanRad,BASeg,GammaFit, gfDelta, gfTau, ...
          NullCondId, SumXtX] = fmri_lddat(datfile)
%
% Load information and parameters from a data file.
% 
% [nNNC,nHEst,DOF,TR,nRuns,nTP,nRows,nCols,nSkip,DTOrder,Rescale,TW,TPS] 
%    = fmri_lddat(datfile) 
%
% datfile - (string) name of data file
%
% nNNC - number of non-null conditions
% nHEst - number of points to estimate in the HDR
% DOF  - degrees of freedom
% TR - time between scans
% nRuns - number of runs
% nTP   - number of scans
% nRows - number of image rows
% nCols - number of image columns
% nSkip - number of initial scans to skip
% DTOrder - order of detrending 
% Rescale - target to which the mean will be rescaled.
% TW - time window
% TPS - prestimulus time
%
% See also fmri_svdat
%
%


%
% fmri_lddat.m
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

if(nargin ~= 1)
  msg = 'Incorrect number of arguments';
  qoe(msg);   error(msg);
end

fid=fopen(deblank(datfile),'r');
if( fid == -1 )
  msg = sprintf('Could not open dof file %s\n',datfile);
  qoe(msg);  error(msg);
end

%% Pre-set selxavg parameters %%
DOF   = 0;
nRuns = 0;
nTP   = 0;
nRows = 0;
nCols = 0;
nSkip = 0;
DTOrder = 0;
Rescale = -1;
HanRad = 0;
BASeg = 0;
nNoiseAC = 0;
GammaFit = 0;
gfDelta = 0;
gfTau = 0;
NullCondId = 0;
SumXtX = 0;

%% --- selavg stuff --- %%
TR      = fscanf2(fid,'%f');
TW      = fscanf2(fid,'%f');
TPS     = fscanf2(fid,'%f');
nPS     = floor(TPS/TR);
nCond   = fscanf2(fid,'%d');
nNNC    = nCond - 1;
nHEst   = fscanf2(fid,'%d'); % nperevent

%% --- selxavg stuff --- %%
DOF   = fscanf2(fid,'%d');
if(isempty(DOF))
  fprintf('INFO: Data file in selavg format');
  return;
end
nRuns = fscanf2(fid,'%d');
nTP   = fscanf2(fid,'%d');
nRows = fscanf2(fid,'%d');
nCols = fscanf2(fid,'%d');
nSkip = fscanf2(fid,'%d');
DTOrder  = fscanf2(fid,'%d');
Rescale  = fscanf2(fid,'%f');
HanRad   = fscanf2(fid,'%f');
nNoiseAC = fscanf2(fid,'%d');
BASeg    = fscanf2(fid,'%d');

%%--------- Gamma Fit -------------%%%
GammaFit   = fscanf2(fid,'%d');
if(isempty(GammaFit))
  GammaFit = 0;
  gfDelta = 0;
  gfTau = 0;
  return;
end
gfDelta  = fscanf2(fid,'%f');
gfTau    = fscanf2(fid,'%f');

NullCondId = fscanf2(fid,'%d');

Nch = nNNC*nHEst;

dummy      = fscanf(fid,'%s',1);
SumXtX     = fscanf(fid,'%f',inf);
nSumXtX = prod(size(SumXtX));
if(nSumXtX ~= Nch*Nch)
  msg = sprintf('Inconsistent Dimensions: nSumXtX=%d, Nch^2 = %d\n',...
                nSumXtX,Nch*Nch);
  qoe(msg);error(msg);
end
fclose(fid);
SumXtX     = reshape(SumXtX,[Nch Nch]);

nPreStim = floor(TPS/TR);


return

%- reads a string followed by a value -----%
function x = fscanf2(fid,fmt)
  fscanf(fid,'%s',1);
  x = fscanf(fid,fmt,1);
return;
