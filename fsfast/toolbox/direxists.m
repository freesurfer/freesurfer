function r = direxists(dirpath)
% r = direxists(dirpath)
%
% returns 1 if dirpath exists, 0 else
%
% $Id: direxists.m,v 1.1 2003/10/02 20:26:57 greve Exp $

r = -1;
if(nargin ~= 1)
  fprintf('r = direxists(dirpath)\n');
  return;
end

d = dir(dirpath);
if(size(d,1) == 0) r = 0;
else               r = 1;
end

return;





