function name = basename(path,ext)
% name = basename(path,<ext>)
% 
% This is an attempt to recreate the unix basename 
% command in matlab. If ext is present, it will check to
% see if ext exists as an extension and strip it if it does.
%
% $Id: basename.m,v 1.4 2011/03/02 00:04:03 nicks Exp $

%
% basename.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:03 $
%    $Revision: 1.4 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

name = [];

if(nargin < 1 | nargin > 2)
  msg = 'name = basename(path,<ext>)'
  return;
end

len = length(path);

if(len == 1) 
  if(strcmp(path,'/'))
    name = '/';
  else
    name = path;
  end
  if(nargin == 2) 
    name = stripextension(name,ext); 
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
    if(nargin == 2) 
      name = stripextension(name,ext); 
    end
    return;
  end
end

name = path;
if(nargin == 2) 
  name = stripextension(name,ext); 
end

return;

%----------------------------------------------------%
function name = stripextension(name,ext)

if(ext(1) ~= '.') 
  ext = ['.' ext]; 
end
next = length(ext);
nname = length(name);
if(nname < next) return; end
tmp = name(nname-next+1:end);
if(~strcmp(tmp,ext)) return; end
name  = name(1:nname-next);
return;
%----------------------------------------------------%



