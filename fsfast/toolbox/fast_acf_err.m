function acferr = fast_acf_err(racf,R,yacf)
% acferr = fast_acf_err(racf,R,yacf)


%
% fast_acf_err.m
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

if(nargin ~= 3)
  fprintf('acferr = fast_acf_err(racf,R,yacf)\n');
  return;
end

nf = length(racf);
Vy = toeplitz(yacf);

% Compute expected racf %
for l = 1:nf
  Dl = diag(ones(nf-l+1,1),l-1);  
  racfexp(l) = trace(R*Dl*R*Vy);
end
racfexp = racfexp/racfexp(1);




