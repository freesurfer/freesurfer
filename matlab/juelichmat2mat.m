function err = juelichmat2mat(matfile)
% err = juelichmat2mat(matfile)
% $Id: juelichmat2mat.m,v 1.1 2007/08/16 02:10:12 greve Exp $

err = 1;
if(nargin ~= 1)
  printf('err = juelichmat2mat(matfile)');
  return;
end

if(~exist(matfile,'file'))
  printf('ERROR: cannot find %s\n',matfile);
  return;
end

bakfile = sprintf('%s.bak',matfile);
if(exist(bakfile,'file'))
  printf('ERROR: backup file already exists for %s\n',matfile);
  return;
end

a = load(matfile);
if(~isfield(a,'M'))
  printf('ERROR: %s does not have an M variable\n',matfile);
  return;
end
M = a.M;

cmd = sprintf('cp %s %s',matfile,bakfile);
s = unix(cmd);
if(s ~= 0)
  printf('ERROR: when executing %s\n',cmd);
  return
end

save(matfile,'M','-v4');

err = 0;
return;

