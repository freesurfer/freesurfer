function mkdirp(path)
% mkdirp(path)
% 
% This is an attempt to recreate the unix "mkdir -p" command
% in matlab.  It has some problems (following links?).
%
% $Id: mkdirp.m,v 1.2 2004/07/16 18:51:32 greve Exp $

if(nargin ~= 1)
  msg = 'USAGE: mkdirp(path)';
  qoe(msg); error(msg);
  return;
end

wd = pwd;

% Get a list of each element in the path %
dirlist = [];
while(~isempty(path) & ~strcmp(path,'/') & ~strcmp(path,'.'))
  base = basename(path);
  dirlist = strvcat(dirlist,base);
  path = dirname(path);
end

% Create the appropriate prefix %
if(strcmp(path,'.'))     pre = './';
elseif(strcmp(path,'/')) pre = '/';
else                     pre = '';
end

% Go backwards (top-bottom) through the path elements 
% and create a directory for each one
for n = size(dirlist,1): -1 : 1
  % dirnm = [dirnm deblank(dirlist(n,:)) '/' ];
  dirnm = [pre deblank(dirlist(n,:))];
  pre = '';

  % Does directory exist? %
  nd = prod(size(dir(dirnm))); 
  if(nd == 0)
    [st msg] = mkdir(dirnm);
    nd = prod(size(dir(dirnm))); 
    if(nd == 0)
      msg = sprintf('ERROR: making dir %s in %s\n',dirnm,pwd);
      if(st==0)  qoe(msg); error(msg); end
    end
  end

  cd(dirnm);
end

cd(wd);

return;

