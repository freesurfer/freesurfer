function name = basename(path)
% name = basename(path)
% 
% This is an attempt to recreate the unix basename 
% command in matlab.
% $Id: basename.m,v 1.1 2003/03/04 20:47:33 greve Exp $

name = [];

if(nargin ~= 1)
  msg = 'USAGE: name = basename(path)'
  qoe(msg); error(msg);
end

len = length(path);

if(len == 1) 
  if(strcmp(path,'/'))
    name = '/';
  else
    name = path;
  end
  return;
end

% strip trailing '/' character %
if(strcmp(path(len),'/'))
  path = path(1:len-1);
  len = len-1;
end

for n = len-1:-1:1
  if(strcmp(path(n),'/'))
    name = path(n+1:len);    
    return;
  end
end

name = path;

return;
