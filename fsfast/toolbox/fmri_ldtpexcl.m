function tTPExclude = fmri_ldtpexcl(TPExclFile)
%
% tTPExclude = ldtpexcl(TPExclFile)
%
% Reads the TP Exclude File.
%
%
%


%
% fmri_ldtpexcl.m
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

tTPExclude = [];

TPExclFile = deblank(TPExclFile);

if( strcmp(TPExclFile,'noexcl')) return; end

fid = fopen(TPExclFile);
if(fid == -1)
  msg = sprintf('Could not open %s',TPExclFile);
  qoe(msg);
  error(msg);
end

tTPExclude = fscanf(fid,'%f');

return;
  
