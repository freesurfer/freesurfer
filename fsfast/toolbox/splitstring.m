function split = splitstring(str)
%
% split = splitstring(str)
%
% Split a string into vertically concatenated strings.
%
% $Id: splitstring.m,v 1.1 2003/03/04 20:47:41 greve Exp $

if(nargin ~= 1)
  msg = 'USAGE: split = splitstring(str)';
  qoe(msg);error(msg);
end

nstr = 1;
[split nscanned ] = sscanf(str,'%s',1);

while(nscanned > 0)
  fmt = repmat('%*s ',[1 nstr]);
  fmt = [fmt '%s'];
  [tmp nscanned] = sscanf(str,fmt,1);
  if(~isempty(tmp))
    split = strvcat(split,tmp);
    nstr = nstr + 1;
  end
end

return
