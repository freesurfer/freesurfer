function hdrdat = fmri_lddat2(datfile)
%
% Load information and parameters from a data file into an hdrdata
% structure. To replace fmri_lddat()
% 
% datfile - (string) name of data file
%
%


%
% fmri_lddat2.m
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
  msg = 'Usage: hdrdat = fmri_lddat2(datfile)'
  qoe(msg); error(msg);
end

fid=fopen(deblank(datfile),'r');
if( fid == -1 )
  msg = sprintf('Could not open dof file %s\n',datfile);
  qoe(msg); error(msg);
end

hdrdat = fmri_hdrdatstruct;

%% --- selavg stuff --- %%
hdrdat.TR      = fscanf2(fid,'%f');
hdrdat.TimeWindow      = fscanf2(fid,'%f');
hdrdat.TPreStim     = fscanf2(fid,'%f');
nPS     = floor(hdrdat.TPreStim/hdrdat.TR);
hdrdat.Nc   = fscanf2(fid,'%d');
hdrdat.Nnnc    = hdrdat.Nc - 1;
hdrdat.Nh   = fscanf2(fid,'%d'); % nperevent

%% --- selxavg stuff --- %%
hdrdat.Version = fscanf2(fid,'%d');

if(isempty(hdrdat.Version))
  %fprintf('    INFO: Data file in selavg format\n');
  hdrdat.Version = 0;
  hdrdat.TER     = hdrdat.TR;
  hdrdat.DOF     = [];
  return;
end

hdrdat.TER   = fscanf2(fid,'%f');
hdrdat.DOF   = fscanf2(fid,'%f');

dummy = fscanf(fid,'%s',1);
hdrdat.Npercond = fscanf(fid,'%d',hdrdat.Nc);

hdrdat.Nruns = fscanf2(fid,'%d');
hdrdat.Ntp   = fscanf2(fid,'%d');
hdrdat.Nrows = fscanf2(fid,'%d');
hdrdat.Ncols = fscanf2(fid,'%d');
hdrdat.Nskip = fscanf2(fid,'%d');
hdrdat.DTOrder  = fscanf2(fid,'%d');
hdrdat.RescaleFactor  = fscanf2(fid,'%f');
hdrdat.HanningRadius   = fscanf2(fid,'%f');
hdrdat.nNoiseAC = fscanf2(fid,'%d');
hdrdat.BrainAirSeg    = fscanf2(fid,'%d');

%%--------- Gamma Fit -------------%%%
hdrdat.GammaFit   = fscanf2(fid,'%d');
if(hdrdat.GammaFit == 0)
  hdrdat.gfDelta = [];
  hdrdat.gfTau = [];
else
  dummy = fscanf(fid,'%s',1);
  hdrdat.gfDelta  = fscanf(fid,'%f',hdrdat.GammaFit);
  dummy = fscanf(fid,'%s',1);
  hdrdat.gfTau    = fscanf(fid,'%f',hdrdat.GammaFit);
end

hdrdat.NullCondId = fscanf2(fid,'%d');

if(~hdrdat.GammaFit)  Nch = hdrdat.Nnnc*hdrdat.Nh;
else                  Nch = hdrdat.Nnnc*length(hdrdat.gfDelta);
end

dummy      = fscanf(fid,'%s',1);
hdrdat.SumXtX     = fscanf(fid,'%f',Nch*Nch);

nSumXtX = prod(size(hdrdat.SumXtX));
if(nSumXtX ~= Nch*Nch )
  msg = sprintf('Inconsistent Dimensions: nSumXtX=%d, Nch^2 = %d\n',...
                nSumXtX,Nch*Nch);
  qoe(msg);error(msg);
end
hdrdat.SumXtX     = reshape(hdrdat.SumXtX,[Nch Nch]);

%--------- Version 2 -------------%
if(hdrdat.Version > 1)
  dummy          = fscanf(fid,'%s',1);
  [hdrdat.hCovMtx c] = fscanf(fid,'%f',Nch*Nch);
  if(c ~= Nch*Nch)
    msg = sprintf('Could not read hCovMtx from %s\n',datfile);
    hdrdat = [];
    qoe(msg); error(msg);
  end
  hdrdat.hCovMtx = reshape(hdrdat.hCovMtx,[Nch Nch]);

  dummy            = fscanf(fid,'%s',1);
  hdrdat.CondIdMap = fscanf(fid,'%f',hdrdat.Nc);
else
  hdrdat.hCovMtx   = hdrdat.SumXtX;
  hdrdat.CondIdMap = [0:hdrdat.Nc];
end

fclose(fid);

return


%%%%%%%%%%%%%%------------------------%%%%%%%%%%%%%%%%%%
%- reads a string followed by a value -----%
function x = fscanf2(fid,fmt)
  fscanf(fid,'%s',1);
  [x c] = fscanf(fid,fmt,1);
  if(c ~= 1)
    x = []; return;
  end
return;
