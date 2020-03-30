function [detvtx,pth] = SStat_mass_FDR2(pval,maskvtx,rate)
% [detvtx,pth] = SStat_mass_FDR2(pval,maskvtx,rate)
%
% Two-stage FDR approach to achieve tighter control of the FDR. This
% procedure is more powerful than the original FDR procedure implemented in
% SStat_mass_FDR.
%
% Input
% pval: P-values.
% maskvtx: Mask's vertices (1-based). Default [] (all vertices included).
% rate: Expected FDR (between 0 and 1).
%
% Output
% detvtx: Detected vertices (1-based).
% pth: FDR threshold.
%
% Original Author: Jorge Luis Bernal Rusiel 
% References: Benjamini, Y., Krieger, A.M., Yekutieli, D. (2006). Adaptive
% linear step-up procedures that control the false discovery rate. 
% Biometrika, 93, 491-507.
%
if nargin < 3
    rate = 0.05;
    if nargin < 2
        maskvtx = [];
    end
end;
nv0 = length(pval);
if isempty(maskvtx)
   maskvtx = 1:nv0; 
end;
p = pval(maskvtx);
nv = length(p);
%% First stage (m0 estimation)
q0 = rate/(1+rate);
pth0 = SStat_mass_FDR(p,q0);
detv0 = maskvtx(p <= pth0);
ndetv0 = length(detv0);
m0 = nv-ndetv0;
%% Second stage
if (ndetv0 ~= 0) && (ndetv0 ~= nv)
    pth = SStat_mass_FDR(p,q0*nv/m0);
    detvtx = maskvtx(p <= pth);
else
    detvtx = detv0;
    pth = pth0;
end;


