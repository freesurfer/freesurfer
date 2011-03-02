function r = fast_fileexists(filename)
% r = fast_fileexists(filename)
% 1 if it exists and is readable , 0 if not


%
% fast_fileexists.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:04 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

if(nargin ~= 1) 
  msg = 'USAGE: r = fast_fileexists(filename)';
  qoe(msg); error(msg);
end

fid = fopen(filename,'r');
if(fid == -1 ) r = 0;
else
  r = 1;
  fclose(fid);
end

return;
