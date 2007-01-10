function C = mkprestimcon(sxa,wcond)
% C = mkprestimcon(sxa,wcond)


%
% mkprestimcon.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
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

C = [];
if(nargin ~= 2)
  fprintf('C = mkprestimcon(sxa,wcond)\n');
  return;
end

nPre = round(sxa.TPreStim/sxa.TER);
nFIR = sxa.Nh;
Nc = sxa.Nc-1;
wcond = wcond(:)';

if(length(wcond) ~= Nc)
  fprintf('ERROR: wcond has wrong number of items\n');
  return;
end

a = zeros(nFIR);
a(:,1:nPre) = -1/nPre;
b = a + eye(nFIR);
C = kron(wcond,b); 

return;
