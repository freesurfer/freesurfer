function slines = readln(fid,nl)
%
% slines = readln(fid,<nlines>)
%
%
% $Id: readln.m,v 1.1 2003/03/04 20:47:41 greve Exp $


F = fread(fid);
lF = length(F);

if( isempty(F) )
  warning('File appears to be empty');
  slines = [];
  return;
end

indNL = find(F == 10); % indicies of new lines
nNL = length(indNL);

if(max(indNL) == lF) % last char is NL
  nLines = nNL;
else
  nLines = nNL + 1;
end

if(nLines == 1)
  slines = setstr(F');
  return;
end

if(nargin == 2) nlmax = nl;
else            nlmax = nLines;
end

if(nlmax > nLines) nlmax = nLines; end

sF = setstr(F');
clear F;
slines = [];
b = 1;
for n = 1:nlmax

  % need this in case file does not end in a newline
  if(n > nNL) e = lF;
  else        e = indNL(n);
  end

  if(e > b + 1)  % skip lines with only NLs
    s = sF(b:e); 
    slines = strvcat(slines,s);
  end

  if(n <= nNL) b = indNL(n) + 1; end

end

return;
