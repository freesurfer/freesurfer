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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:06 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
  
