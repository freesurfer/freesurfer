function [curv, fnum] = fast_read_curv(fname)
%
% [curv, fnum] = fast_read_curv(fname)
% reads a binary curvature file into a vector
%

% open it as a big-endian file
fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
 str = sprintf('could not open curvature file %s.', fname) ;
 error(str) ;
end

vnum = fread3(fid) ;
NEW_VERSION_MAGIC_NUMBER = 16777215;
if (vnum == NEW_VERSION_MAGIC_NUMBER)
  vnum = fread(fid, 1, 'int32') ;
  fnum = fread(fid, 1, 'int32') ;
  vals_per_vertex = fread(fid, 1, 'int32') ;
  curv = fread(fid, vnum, 'float') ; 
  fclose(fid) ;
else
  fnum = fast_fread3(fid) ;
  curv = fread(fid, vnum, 'int16') ./ 100 ; 
  fclose(fid) ;
end


