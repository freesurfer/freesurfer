function [item, count] = sscanfitem(tline,nthitem)
% [item count] = sscanfitem(tline,nthitem)
%
% Reads the nth item from a string list of items separated by white
% space. The item is returned as a string. If successful, count=1,
% otherwise count=0 (and item is empty).
%
%


%
% sscanfitem.m
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

item = '';
count = 0;

if(nargin ~= 2)
  fprintf('[item count] = sscanfitem(tline,nthitem)\n');
  return;
end

fmt = '%s';
for n = 1:nthitem-1
  fmt = sprintf('%%*s %s',fmt);
end

[item count] = sscanf(tline,fmt,1);

return;





