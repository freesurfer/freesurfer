function [FPR, alpha] = ComputeFPR(p, pDelta)
%
% [FPR alpha] = ComputeFPR(p, <pDelta>)
%

if(nargin == 1)  pDelta = .01; end

pRange = [(min(p)+pDelta/2) pDelta  (max(p)-pDelta/2)];
[pDist0 alpha pDelta0] = samp2pdf(p,pRange);

FPR = cumtrapz(alpha,pDist0);
FPR = FPR + alpha(1);

return;
