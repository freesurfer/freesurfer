function w = fast_read_wfile(fname)
%
% w = fast_read_wfile(fname)
% reads a vector from a binary 'w' file
%	fname - name of file to read from
%	w     - vector of values 
%


%
% fast_read_wfile.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
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

w = [];
v = [];

if(nargin ~= 1)
  fprintf('USAGE: [w,v] = fast_read_wfile(fname, w) \n');
  return;
end

% open it as a big-endian file
fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
  str = sprintf('could not open w file %s.', fname) ;
  error(str) ;
end

fread(fid, 1, 'int16') ;  % Skip ilat
vnum = fast_fread3(fid) ; % Number of non-zero values
v = zeros(vnum,1) ;
w0 = zeros(vnum,1) ;
for i=1:vnum
  v(i) = fast_fread3(fid) ;
  w0(i) = fread(fid, 1, 'float') ;
end

fclose(fid) ;

w = zeros(max(v),1);
w(v+1) = w0;

return;




