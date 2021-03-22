function runlist = fast_runlistfile(runlistfile)
% runlist = fast_runlistfile(runlistfile)


%
% fast_runlistfile.m
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

runlist = [];

if(nargin ~= 1)
  msg = 'runlist = fast_runlistfile(runlistfile)';
  qoe(msg);error(msg);
end

fid = fopen(runlistfile);
if(fid == -1)
  msg = sprintf('Could not open %s',runlistfile);
  qoe(msg);error(msg);
end

runid = deblank(fscanf(fid,'%s',1));
while( ~isempty(runid) )
  runlist = strvcat(runlist,runid);
  runid = fscanf(fid,'%s',1);
end
fclose(fid);

if(isempty(runlist))
  msg = sprintf('ERROR: no runs found in %s\n',runlistfile);
  qoe(msg);error(msg);
end

return;
