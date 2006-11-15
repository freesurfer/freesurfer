function st = fast_ldstf(stf)
% st = fast_ldstf(stf)
%
% Load a stimulus timing file - 2 or 3 cols
% Col 1: stimulus onset time (sec)
% Col 2: stimulus duration (sec)
% Col 3: stimulus weight
%
% Quitely returns an empty if file does not exist.
%
% $Id: fast_ldstf.m,v 1.3 2006/11/15 22:50:44 greve Exp $

st = [];

if(nargin ~= 1)
  fprintf('st = fast_ldstf(stf)\n');
  return;
end

fp = fopen(stf,'r');
if(fp == -1)
  %fprintf('ERROR: could not open %s\n',stf);
  % Don't print an error here, in case allowing some conditions to 
  % not exist in some runs. Just return an empty.
  return;
end

use3col = 0;
nthpres = 1;  
while(1)

  % scroll through any blank lines or comments
  while(1)
    tline = fgetl(fp);
    if(~isempty(tline) & tline(1) ~= '#') break; end
  end
  if(tline(1) == -1) break; end

  [item c] = sscanfitem(tline,1);
  if(c ~= 1) fprintf('Format error %s\n',stf); st=[]; return; end
  st(nthpres,1) = sscanf(item,'%f',1); 
  
  [item c] = sscanfitem(tline,2);
  if(c ~= 1) fprintf('Format error %s\n',stf); st=[]; return; end
  st(nthpres,2) = sscanf(item,'%f',1); 
  
  [item c] = sscanfitem(tline,3);
  if(c == 1) 
    st(nthpres,3) = sscanf(item,'%f',1); 
    use3col = 1;
  else
    if(use3col)
      fprintf('Format error %s (3 col)\n',stf);
      st=[];
      return;
    end
  end

  nthpres = nthpres + 1;  

end % while (1)

fclose(fp);


return;
