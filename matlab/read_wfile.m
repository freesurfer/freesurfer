function [w,v] = read_wfile(fname)
%
% [w,v] = read_wfile(fname)
% reads a vector into a binary 'w' file
% fname - name of file to write to
% w     - vector of values to be written
% v     - vector of vertex indices (0-based)
%
% See also write_wfile.
%
% $Id: read_wfile.m,v 1.2 2003/08/20 17:22:13 greve Exp $

w = [];
v = [];

if(nargin ~= 1)
  fprintf('[w,v] = read_wfile(fname)\n');
  return;
end

% open it as a big-endian file
fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
  str = sprintf('could not open w file %s.', fname) ;
  error(str) ;
end

fread(fid, 1, 'int16') ;
vnum = fread3(fid) ;
w = zeros(vnum,1) ;
v = zeros(vnum,1) ;
for i=1:vnum
  v(i) = fread3(fid) ; % vertex number (0-based)
  w(i) = fread(fid, 1, 'float') ; % vertex value
end

fclose(fid) ;






