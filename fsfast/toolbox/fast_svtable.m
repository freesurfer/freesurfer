function err = fast_svtable(tbl, rowid, colid, fname, base)
% err = fast_svtable(tbl, rowid, colid, fname, <base>)
%
% Prints out a data table that is compatible with fast_ldtable.
% The first row is a list of column headers (colid), the first
% column is a list of row headers (rowid). The first field (ie,
% string in the first row and col) will be base (if no base is
% supplied then 0000 is put as a place-holder).
%
% $Id: fast_svtable.m,v 1.1 2003/05/28 20:55:59 greve Exp $

err = 1;

if(nargin < 4 | nargin > 5)
  fprintf('err = fast_svtable(tbl, rowid, colid, fname, <base>)\n');
  return;
end

if(exist('base') ~= 1) base = []; end
if(isempty(base)) base = '0000'; end

sztbl = size(tbl);
if(sztbl(1) ~= size(rowid,1))
  fprintf('Table size row dimension mismatch\n');
  return;
end
if(sztbl(2) ~= size(colid,1))
  fprintf('Table size col dimension mismatch\n');
  return;
end

fid = fopen(fname,'w');
if(fid == -1)
  fprintf('ERROR: could not open %s for writing\n',fname);
  return;
end

for r = 0:size(rowid,1);
  for c = 0:size(colid,1)
    if(r==0 & c ==0) 
      fprintf(fid,'%s',base);
      continue;
    end
    if(r==0) % note: c cannot = 0 here
      fprintf(fid,' %s',deblank(colid(c,:)));
      continue;
    end
    if(c==0)  % note: r cannot = 0 here
      fprintf(fid,'%s',deblank(rowid(r,:)));
      continue;
    end
    fprintf(fid,' %g',tbl(r,c));
  end
  fprintf(fid,'\n');
end
fprintf(fid,'\n');

err = 0;

return;
