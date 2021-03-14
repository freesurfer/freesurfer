function [FPR, alpha] = ComputeFPR(p, pDelta)
%
% [FPR alpha] = ComputeFPR(p, <pDelta>)
%


%
% ComputeFPR.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

if(nargin == 1)  pDelta = .01; end

pRange = [(min(p)+pDelta/2) pDelta  (max(p)-pDelta/2)];
[pDist0 alpha pDelta0] = samp2pdf(p,pRange);

FPR = cumtrapz(alpha,pDist0);
FPR = FPR + alpha(1);

return;
