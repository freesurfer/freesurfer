function corslice = fast_ldcorslice(corslicefile);
%
% corslice = fast_ldcorslice(corslicefile);
%


%
% fast_ldcorslice.m
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
