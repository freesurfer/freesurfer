function [tbl, rowid, colid] = fast_ldtable(tablefile)
% [tbl, rowid, colid] = fast_ldtable(tablefile)
%
% Loads a text data table where the data table is 
% formatted in the following way:
%   1. The first row is a list of strings that identify
%      each column (returned in colid)
%   2. The first col is a list of strings that identify
%      echo row (returned in rowid)
%   3. The first item is ignored (but must be there)
%   4. Blank lines are ignored
%   5. Any line begining with a # is ignored
%
% $Id: fast_ldtable.m,v 1.2 2004/11/14 22:32:46 greve Exp $
%

tbl=[];
colid=[];
rowid=[];

if(nargin ~= 1)
  fprintf('[tbl, rowid, colid] = fast_ldtable(tablefile)\n');
  return;
end

fid = fopen(tablefile);
if(fid == -1)
  fprintf('ERROR: opening %s\n',tablefile);
  return;
end

tline = fgetl(fid);
if(tline == -1)
  fprintf('ERROR: %s is not correctly formatted, no first line\n', ...
	  tablefile);
  fclose(fid);
  return;
end
colid = splitstring(tline);
colid = colid(2:end,:);

%----------- Loop through all the lines ----------------------%
nthrow = 1;
while(1)

  % scroll through any blank lines or comments %
  while(1)
    tline = fgetl(fid);
    if(~isempty(tline) & tline(1) ~= '#') break; end
  end
  if(tline(1) == -1) break; end

  % Get count of number of items in the line %
  [items count] = sscanf(tline,'%s');
  if(nthrow == 1) count0 = count;
  else
    if(count == count0)
      fprintf('ERROR: %s is not correctly formatted.\n',tablefile);
      fprintf(' Line %d has %d items, expecting %d\n',...
	      nthrow,count,count0);
      fclose(fid);
      return;
    end
  end
  
  % Read the row id %  
  key = sscanf(tline,'%s',1);
  rowid = strvcat(rowid,key);

  % Read the values in the row %
  rowvals = [];
  for c = 1:count
    valcol = str2num(sscanfitem(tline,c+1));
    rowvals = [rowvals valcol];
  end
  tbl = [tbl; rowvals];
end % while (1)

fclose(fid);

return;
%---------------------------------------------------------%
