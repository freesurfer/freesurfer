function racfexp = fast_racf_exp(R,yacf)
% racfexp = fast_racf_exp(R,yacf)
%
% Computes the expected value of the acf of the residual
% given the acf of y and the residual error forming matrix R.
%
% Based on Worsley, et al, NeuroImage 15, 1-15, 2002.
%


%
% fast_racf_exp.m
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

if(nargin ~= 2)
  fprintf('racfexp = fast_racf_exp(R,yacf)\n');
  return;
end

nf = length(yacf);
Vy = toeplitz(yacf);
RVy = R*Vy;

for l = 1:nf
  Dl = diag(ones(nf-l+1,1),l-1);  
  racfexp(l) = trace(R*Dl*RVy);
end
racfexp = racfexp/racfexp(1);

return;


