function corslice = fast_ldcorslice(corslicefile);
%
% corslice = fast_ldcorslice(corslicefile);
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
