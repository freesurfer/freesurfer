function [FPR, alpha] = ComputeFPR(p, pDelta)
%
% [FPR alpha] = ComputeFPR(p, <pDelta>)
%


%
% ComputeFPR.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:29 $
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

if(nargin == 1)  pDelta = .01; end

pRange = [(min(p)+pDelta/2) pDelta  (max(p)-pDelta/2)];
[pDist0 alpha pDelta0] = samp2pdf(p,pRange);

FPR = cumtrapz(alpha,pDist0);
FPR = FPR + alpha(1);

return;
