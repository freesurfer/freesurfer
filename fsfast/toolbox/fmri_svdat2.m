function fmri_svdat(datfile,hd,format);
%
% fmri_svdat(datfile,hdrdat,<format>)
%
%


%
% fmri_svdat2.m
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

if(nargin ~= 2 & nargin ~= 3)
  msg = 'Usage: fmri_svdat(datfile,hdrdat,<format>)';
  qoe(msg); error(msg);
end

if(nargin == 2) format = 'selxavg'; end

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

if(strcmp(lower(format),'selavg')) return; end

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
fprintf(fid,'%g\n',hd.SumXtX);

%-------- Version 2 ------------- %
if(hd.Version > 1)

  fprintf(fid,'hCovMtx\n');
  fprintf(fid,'%g\n',hd.hCovMtx);

  fprintf(fid,'CondIdMap ');
  fprintf(fid,'%3d ',hd.CondIdMap);  
  fprintf(fid,'\n');

end


fclose(fid);



