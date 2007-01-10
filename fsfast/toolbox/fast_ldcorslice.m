function corslice = fast_ldcorslice(corslicefile);
%
% corslice = fast_ldcorslice(corslicefile);
%


%
% fast_ldcorslice.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

if(nargin ~= 1)
  msg = 'USAGE: corslice = fast_ldcorslice(corslicefile);';
  qoe(msg);error(msg);
end

%%%% Open the corslicefile %%%%%
Endian = 0;
if(Endian == 0) fid=fopen(corslicefile,'r','b'); % Big-Endian
else            fid=fopen(corslicefile,'r','l'); % Little-Endian
end
if(fid == -1)
  msg = sprintf('Could not open %s for reading.',corslicefile); 
  qoe(msg); error(msg);
end

%%% Read the file in corslicefile %%%
precision = 'uint8';
Nv = 256*256;
z = fread(fid,Nv,precision);
corslice = reshape(z, [256 256])'; %' transpose for row major
fclose(fid); 


return;
