function [item, count] = sscanfitem(tline,nthitem)
% [item count] = sscanfitem(tline,nthitem)
%
% Reads the nth item from a string list of items separated by white
% space. The item is returned as a string. If successful, count=1,
% otherwise count=0 (and item is empty).
%
% $Id: sscanfitem.m,v 1.1 2004/10/16 05:27:25 greve Exp $

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





