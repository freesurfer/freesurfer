function fmri_svsfa(sfa,stem)
% fmri_svsfa(sfa,stem)


%
% fmri_svsfa.m
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
  msg = 'USAGE: fmri_svsfa(sfa,stem)';
  qoe(msg); error(msg);
end

stem = deblank(stem);

%% Open the sfa file %%
sfaname = sprintf('%s.sfa',stem);
fid = fopen(sfaname,'w');
if(fid == -1) 
  msg = sprintf('Could not open %s for writing',sfaname);
  qoe(msg); error(msg);
end

fprintf(fid,'SelectiveFrequencyAverage\n');
svstruct(sfa,fid);
fclose(fid);

return;




fprintf(fid,'version     %d\n',sfa.version);
fprintf(fid,'type        %s\n',sfa.type);
fprintf(fid,'dof         %d\n',sfa.dof);
fprintf(fid,'TR          %f\n',sfa.TR);
fprintf(fid,'Ntp         %d\n',sfa.Ntp);
fprintf(fid,'ncycles     %d\n',sfa.ncycles);
fprintf(fid,'fundamental %f\n',sfa.fundamental);
fprintf(fid,'delay       %f\n',sfa.delay);
fprintf(fid,'delay_stem  %s\n',sfa.delay_stem);
fprintf(fid,'direction   %s\n',sfa.direction);
fprintf(fid,'freqskip    %d\n',sfa.freqskip);
fprintf(fid,'skirtskip   %d\n',sfa.skirtskip);

fprintf(fid,'nsignal %d\n',length(sfa.isignal));
fprintf(fid,'isignal ');  
fprintf(fid,' %3d',sfa.isignal);  
fprintf(fid,'\n');

fprintf(fid,'nnoise  %d\n',length(sfa.inoise));
fprintf(fid,'inoise ');   
fprintf(fid,' %3d',sfa.inoise);   
fprintf(fid,'\n');

fprintf(fid,'nexclude %d\n',length(sfa.iexclude));
fprintf(fid,'iexclude '); 
fprintf(fid,' %3d',sfa.iexclude); 
fprintf(fid,'\n');

fprintf(fid,'nrows %d\n',sfa.nrows);
fprintf(fid,'ncols %d\n',sfa.ncols);

fprintf(fid,'nslices %d\n',length(sfa.slice_delay));
fprintf(fid,'slice_delay ');
fprintf(fid,'%f ',sfa.slice_delay);
fprintf(fid,'\n');

fprintf(fid,'meanval %f\n',sfa.meanval);
fprintf(fid,'rescale_target %f\n',sfa.rescale_target);
fprintf(fid,'rescale_factor %f\n',sfa.rescale_factor);
fprintf(fid,'nksip %d\n',sfa.nskip);
fprintf(fid,'hanrad %f\n',sfa.hanrad);
fprintf(fid,'detrend %d\n',sfa.detrend);

fprintf(fid,'nactiveharm %d\n',length(sfa.slice_delay));
fprintf(fid,'nactiveharm ');
fprintf(fid,'%d ',sfa.activeharm);
fprintf(fid,'\n');


fclose(fid);

return
