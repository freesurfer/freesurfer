function tTPExclude = fmri_ldtpexcl(TPExclFile)
%
% tTPExclude = ldtpexcl(TPExclFile)
%
% Reads the TP Exclude File.
%
% $Id: fmri_ldtpexcl.m,v 1.1 2003/03/04 20:47:40 greve Exp $
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
  
