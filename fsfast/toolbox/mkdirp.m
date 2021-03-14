function err = mkdirp(path)
% err = mkdirp(path)
% 
% This is an attempt to recreate the unix "mkdir -p" command
% in matlab.  It has some problems (following links?).
% Now I just use a unix() command
%
%


%
% mkdirp.m
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

err = 1;

if(nargin ~= 1)
  msg = 'USAGE: mkdirp(path)';
  qoe(msg); error(msg);
  return;
end

cmd = sprintf('mkdir -p %s',path);
err = unix(cmd);
if(err) fprintf('ERROR: %s\n',cmd); end
return;

%---------------------------------------------------------%
% Never gets down here


wd = pwd;

% Get a list of each element in the path %
dirlist = [];
while(~isempty(path) & ~strcmp(path,'/') & ~strcmp(path,'.'))
  base = basename(path);
  dirlist = strvcat(dirlist,base);
  path = fast_dirname(path);
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

