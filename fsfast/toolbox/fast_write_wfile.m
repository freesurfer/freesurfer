function [w] = fast_write_wfile(fname, w)
%
% [w] = fast_write_wfile(fname, w)
% writes a vector into a binary 'w' file
%  fname - name of file to write to
%  w     - vector of values to be written
%

if(nargin ~= 2)
  fprintf('USAGE: [w] = fast_write_wfile(fname, w) \n');
  return;
end

% open it as a big-endian file
fid = fopen(fname, 'wb', 'b') ;
vnum = length(w) ;

fwrite(fid, 0, 'int16') ;
fast_fwrite3(fid, vnum) ;
for i=1:vnum
  fast_fwrite3(fid, i-1) ;
  wt = w(i) ;
  fwrite(fid, wt, 'float') ;
end

fclose(fid) ;

return
