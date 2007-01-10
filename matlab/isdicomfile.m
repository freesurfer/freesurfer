function [res, dcminfo] = isdicomfile(fname)
% [res, dcminfo] = isdicomfile(fname)
%
% Determines whether the given file is dicom by
% trying to read the dicom info. If this fails,
% then res=0 and dcminfo=[]. If successful, then
% res=1, and dcminfo is the result of matlabs
% dicominfo.
%


%
% isdicomfile.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
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
