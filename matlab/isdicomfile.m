function [res, dcminfo] = isdicomfile(fname)
% [res, dcminfo] = isdicomfile(fname)
%
% Determines whether the given file is dicom by
% trying to read the dicom info. If this fails,
% then res=0 and dcminfo=[]. If successful, then
% res=1, and dcminfo is the result of matlabs
% dicominfo.
%
% $Id: isdicomfile.m,v 1.1 2003/01/23 20:01:34 greve Exp $

try
  dcminfo = dicominfo(fname);
  res = 1;
catch
  dcminfo = [];
  res = 0;
end

return;

%%%%%%% This was Anders original code, does not always work %%%%%%%%
fid = fopen(fname,'r');
if fid < 0
  res = 0;
else
  stat = fseek(fid,128,'bof'); % move to DICM string
  tmp = char(fread(fid,4,'uchar')');%'
  res = strcmp(tmp,'DICM');
  fclose(fid);
end
