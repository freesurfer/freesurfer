function val = tdr_measasc(measasc,varname)
% val = tdr_measasc(measasc,varname)
%
% $Id: tdr_measasc.m,v 1.1 2003/11/06 19:46:27 greve Exp $

val = [];
if(nargin ~= 2)
  fprintf('val = tdr_measasc(measasc,varname)\n');
  return;
end

fp = fopen(measasc,'r');
if(fp == -1)
  fprintf('ERROR: could not open %s\n',measasc);
  return;
end

while(1)

  % scroll through any blank lines or comments %
  while(1)
    tline = fgetl(fp);
    if(~isempty(tline) & tline(1) ~= '#') break; end
  end
  if(tline(1) == -1) break; end

  % Get count of number of items in the line %
  [items count] = sscanf(tline,'%s');

  % Read the keyword %  
  key = sscanf(tline,'%s',1);
  %fprintf('key = %s\n',key);

  if(strcmp(key,varname))
    valstr = sscanf(tline,'%*s %*s %s',1);
    ind = find(valstr ~= '[' & valstr ~= ']');
    valstr = valstr(ind);
    val = sscanf(valstr,'%f');
    %keyboard
    return;
  end

end % while (1)

fprintf('ERROR: could not fine %s in %s\n',varname,measasc);


return;