function fast_fwrite3(fid, val)
% fast_fwrite3(fid, val)
%
% Write a 3 byte integer into a file
%


%
% fast_fwrite3.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
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
  fprintf('USAGE: fast_fwrite3(fid, val)\n');
  return;
end

%fwrite(fid, val, '3*uchar') ;
b1 = bitand(bitshift(val, -16), 255) ;
b2 = bitand(bitshift(val, -8), 255) ;
b3 = bitand(val, 255) ; 
fwrite(fid, b1, 'uchar') ;
fwrite(fid, b2, 'uchar') ;
fwrite(fid, b3, 'uchar') ;

return
