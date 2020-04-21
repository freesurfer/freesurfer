function retpar = fast_ldretpar(retparfile)
% retpar = fast_ldretpar(retparfile)
%
% Reads in a retinotopy paradigm file
% retpar.direction = direction;
% retpar.stimtype = stimtype;
% 
% File should look something like:
%   stimtype eccen
%   direction pos
%

retpar = [];
if(nargin ~= 1)
  fprintf('retpar = fast_ldretpar(retparfile)\n');
  return;
end

if(~exist(retparfile,'file'))
  fprintf('ERROR: cannot find %s\n',retparfile);
  return;
end

stimtype  = '';
direction = '';
fp = fopen(retparfile,'r');
if(fp == -1)
  fprintf('ERROR: could not open %s\n',retparfile);
  return;
end
nthline = 1;
while(1)
  % scroll through any blank lines or comments
  while(1)
    tline = fgetl(fp);
    if(~isempty(tline) & tline(1) ~= '#') break; end
  end
  if(tline(1) == -1) break; end

  key = sscanf(tline,'%s',1);
  %fprintf('key = %s\n',key);
  
  switch(key)
   case 'stimtype',    stimtype  = sscanf(tline,'%*s %s',1);
   case 'direction',   direction = sscanf(tline,'%*s %s',1);
   otherwise
    fprintf('INFO: key %s unrecognized, line %d, skipping\n',key,nthline);
  end
  nthline = nthline + 1;
end % while (1)
fclose(fp);

if(isempty(stimtype))
  fprintf('ERROR: retinotopy paradigm file %s not formatted correctly\n',retparfile);
  fprintf('It is missing a stimtype field\n');
  return;
end
if(isempty(direction))
  fprintf('ERROR: retinotopy paradigm file %s not formatted correctly\n',retparfile);
  fprintf('It is missing a direction field\n');
  return;
end

if(~strcmp(stimtype,'eccen') & ~strcmp(stimtype,'polar'))
  fprintf('ERROR: retinotopy paradigm file %s not formatted correctly\n',retparfile);
  fprintf('The stimtype field is set to %s, must be eccen or polar\n',stimtype);
  return;
end  
if(~strcmp(direction,'pos') & ~strcmp(direction,'neg'))
  fprintf('ERROR: retinotopy paradigm file %s not formatted correctly\n',retparfile);
  fprintf('The direction field is set to %s, must be pos or neg\n',direction);
  return;
end  

retpar.direction = direction;
retpar.stimtype = stimtype;

return;

