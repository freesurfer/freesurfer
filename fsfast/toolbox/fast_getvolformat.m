function volformat = fast_getvolformat(volid)
% volformat = fast_getvolformat(volid)
% 
% $Id: fast_getvolformat.m,v 1.1 2003/03/04 20:47:38 greve Exp $

volformat = [];

if(nargin ~= 1)
  msg = 'USAGE: volformat = fast_getvolformat(volid)'
  qoe(msg); error(msg);
end

if(fast_ismincvol(volid))  volformat = 'minc';
elseif(fast_isbvol(volid)) volformat = 'bfile';
end


return;
