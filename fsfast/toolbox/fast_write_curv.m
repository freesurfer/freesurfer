function [curv] = write_curv(fname, curv, fnum)
%
% [curv] = write_curv(fname, curv, fnum)
%
% writes a curvature vector into a binary file
% fname - name of file to write to
% curv  - vector of curvatures
% fnum  - # of faces in surface.
%

if(nargin ~= 3)
  fprintf('USAGE: curv = write_curv(fname, curv, fnum)\n');
  return;
end

% open it as a big-endian file
fid = fopen(fname, 'wb', 'b') ;
vnum = length(curv) ;
NEW_VERSION_MAGIC_NUMBER = 16777215;
fast_fwrite3(fid, NEW_VERSION_MAGIC_NUMBER ) ;
fwrite(fid, vnum,'int32') ;
fwrite(fid, fnum,'int32') ;
fwrite(fid, 1, 'int32');
fwrite(fid, curv, 'float') ;
fclose(fid) ;

