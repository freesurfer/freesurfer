function z = fast_p2z(p)
% z = fast_p2z(p)
%
% Converts a significance to a z-score. p must be
% between -1 and 1
%
%
%


%
% fast_p2z.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:04 $
%    $Revision: 1.4 $
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

z = [];
if(nargin ~= 1) 
  fprintf('z = fast_p2z(p)\n');
  return;
end

z = sqrt(2)*erfcinv(2*abs(p)) .* sign(p);


return;










