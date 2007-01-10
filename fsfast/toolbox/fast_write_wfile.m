function [w] = fast_write_wfile(fname, w)
%
% [w] = fast_write_wfile(fname, w)
% writes a vector into a binary 'w' file
%  fname - name of file to write to
%  w     - vector of values to be written
%


%
% fast_write_wfile.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
%    $Revision: 1.2 $
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
