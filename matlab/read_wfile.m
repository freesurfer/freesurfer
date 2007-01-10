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


%
% read_wfile.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:10 $
%    $Revision: 1.3 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%


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






