function fmri_svdat3(datfile,hd)
%
% fmri_svdat3(datfile,hdrdat)
%
%


%
% fmri_svdat3.m
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

if(nargin ~= 2)
  msg = 'Usage: fmri_svdat(datfile,hdrdat)';
  qoe(msg); error(msg);
end

if(hd.Version == 2)
  fmri_svdat2(datfile,hd);
  return;
end

if(isempty(hd.runlist)) 
   hd.runlist = -ones(hd.Nruns,1);
else
  if(length(hd.runlist) ~= hd.Nruns)
    msg = sprintf('ERROR: number of runs does not equal run list\n');
    qoe(msg);  error(msg);
  end
end
if(isempty(hd.funcstem))  hd.funcstem = 'unknown'; end
if(isempty(hd.parname))   hd.parname  = 'unknown'; end
if(isempty(hd.extregstem))
   hd.extregstem = 'none';
   hd.nextreg = 0;
   hd.extregortho = 0;
end
if(isempty(hd.nextreg)) hd.nextreg = -1; end
if(isempty(hd.extregortho)) hd.extregortho = 0; end

%% Open the output data file %%
fid=fopen(deblank(datfile),'w');
if( fid == -1 )
  msg = sprintf('Could not open dof file %s\n',datfile);
  qoe(msg);  error(msg);
end

%% --- selavg stuff --- %%
fprintf(fid,'TR         %g\n',hd.TR);
fprintf(fid,'TimeWindow %g\n',hd.TimeWindow);
fprintf(fid,'TPreStim   %g\n',hd.TPreStim);
fprintf(fid,'nCond      %d\n',hd.Nc); % includes fixation
fprintf(fid,'Nh        %3d\n',hd.Nh);

%% --- selxavg stuff --- %%
fprintf(fid,'Version    %g\n',hd.Version);
fprintf(fid,'TER        %g\n',hd.TER);

if(~hd.GammaFit) hd.Nh  = floor(hd.TimeWindow/hd.TR);
else             hd.Nh  = 1;
end

Nch = hd.Nh*hd.Nnnc;
fprintf(fid,'DOF      %3d\n',hd.DOF);
fprintf(fid,'Npercond ');
fprintf(fid,'%d ',hd.Npercond);
fprintf(fid,'\n');
fprintf(fid,'nRuns   %3d\n',hd.Nruns);
fprintf(fid,'nTP     %3d\n',hd.Ntp);
fprintf(fid,'Rows    %3d\n',hd.Nrows);
fprintf(fid,'Cols    %3d\n',hd.Ncols);
fprintf(fid,'nSkip   %3d\n',hd.Nskip);
fprintf(fid,'DTOrder %3d\n',hd.DTOrder);
fprintf(fid,'Rescale  %g\n',hd.RescaleFactor);
fprintf(fid,'HanRad   %g\n',hd.HanningRadius);
fprintf(fid,'nNoiseAC %3d\n',hd.nNoiseAC);
fprintf(fid,'BASeg    %d\n',hd.BrainAirSeg);

fprintf(fid,'GammaFit %d\n',hd.GammaFit);
if(hd.GammaFit > 0) 
  fprintf(fid,'gfDelta ');
  fprintf(fid,'%g ',hd.gfDelta);
  fprintf(fid,'\n');
  fprintf(fid,'gfTau ');
  fprintf(fid,'%g ',hd.gfTau);
  fprintf(fid,'\n');
end

fprintf(fid,'NullCondId %d\n',hd.NullCondId);
fprintf(fid,'SumXtX\n');
if(isempty(hd.SumXtX)) hd.SumXtX = eye(Nch); end
fprintf(fid,'%g\n',hd.SumXtX);

fprintf(fid,'hCovMtx\n');
if(isempty(hd.hCovMtx)) hd.hCovMtx = eye(Nch); end
fprintf(fid,'%g\n',hd.hCovMtx);

fprintf(fid,'CondIdMap ');
fprintf(fid,'%3d ',hd.CondIdMap);  
fprintf(fid,'\n');

%---------- Version 3 stuff ---------------------%
fprintf(fid,'LPFFlag %d\n',hd.LPFFlag);
fprintf(fid,'HPF     %d %d\n',hd.HPF(1),hd.HPF(2));
fprintf(fid,'WhitenFlag %d\n',hd.WhitenFlag);
fprintf(fid,'RunList ');
fprintf(fid,'%3d ',hd.runlist);
fprintf(fid,'\n');
fprintf(fid,'FuncStem %s \n',hd.funcstem);
fprintf(fid,'ParName  %s \n',hd.parname);
fprintf(fid,'ExtRegStem  %s \n',hd.extregstem);
fprintf(fid,'NExtReg %d \n',hd.nextreg);
fprintf(fid,'ExtRegOrtho %d \n',hd.extregortho);

fclose(fid);



