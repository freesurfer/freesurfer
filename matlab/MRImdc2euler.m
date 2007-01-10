function eulerangles = MRImdc2euler(Mdc)
% eulerangles = MRImdc2euler(Mdc)
%
% Solves for the Euler angles given the 3x3 matrix of direction cosines.
% This code does not work in all cases.
%


%
% MRImdc2euler.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
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


eulerangles = [];
if(nargin ~= 1)
  fprintf('eulerangles = MRImdc2euler(Mdc)\n');
  return;
end

theta = acos(Mdc(3,3));
phi = asin(Mdc(3,2)/(sin(theta)+eps));
psi = asin(Mdc(2,3)/(sin(theta)+eps));

eulerangles = [theta phi psi];

return;



