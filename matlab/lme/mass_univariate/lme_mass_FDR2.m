function [detvtx,sided_pval,pth,m0] = lme_mass_FDR2(pval,sgn,maskvtx,rate,tail)
% [detvtx,sided_pval,pth,m0] = lme_mass_FDR2(pval,maskvtx,rate,tail)
%
% Two-stage FDR approach to achieve tighter control of the FDR. This
% procedure is more powerful than the original FDR procedure implemented in
% lme_mass_FDR.
%
% Input
% pval: P-values.
% sgn: Sign of vertex-wise contrasts.
% maskvtx: Mask's vertices (1-based). Default [] (all vertices included).
% rate: Expected FDR.
% tail: Must be -1 for left-sided or 0 for two-sided or 1 for right-sided 
% hypothesis testing.
%
% Output
% detvtx: Detected vertices (1-based).
% spval: Unsigned sided p-values.
% pth: FDR threshold.
% m0: Estimated number of null vertices.
%
% $Revision: 1.1.2.2 $  $Date: 2013/02/23 21:08:10 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/02/23 21:08:10 $
%    $Revision: 1.1.2.2 $
% References: Benjamini, Y., Krieger, A.M., Yekutieli, D. (2006). Adaptive
% linear step-up procedures that control the false discovery rate. 
% Biometrika, 93, 491-507.
%
if nargin < 2
    error('Too few inputs');
elseif nargin < 5
    tail = 0;
    if nargin < 4
        rate = 0.05;
        if nargin < 3
            maskvtx = [];
        end
    end;
end;
nv0 = length(pval);
if isempty(maskvtx)
   maskvtx = 1:nv0; 
end;
p = pval(maskvtx).*sgn(maskvtx);
p(pval(maskvtx)==1) = 1;
nv = length(p);

%% First stage (m0 estimation)
q0 = rate/(1+rate);
spval = sidedPval(p,tail);
pth0 = lme_mass_FDR(spval,q0);
if tail==0
    % two-sided thresh
     detv0 = maskvtx(spval <= pth0);
elseif  tail==-1
    % left-sided thresh
    vtx = maskvtx(p<0);
    detv0 = vtx(spval(p<0) <= pth0);
elseif tail==1
    % right-sided thresh
    vtx = maskvtx(p>0);
    detv0 = vtx(spval(p>0) <= pth0);
end;
ndetv0 = length(detv0);
m0 = nv-ndetv0;
%% Second stage
if (ndetv0 ~= 0) && (ndetv0 ~= nv)
    % one sided-thresh
    pth = lme_mass_FDR(spval,q0*nv/m0);
    if tail==0
        % two-sided thresh
        detvtx = maskvtx(spval <= pth);
    elseif  tail==-1
        % left-sided thresh
        detvtx = vtx(spval(p<0) <= pth);
    elseif tail==1
        % right-sided thresh
        detvtx = vtx(spval(p>0) <= pth);
    end;
else
    detvtx = detv0;
    pth = pth0;
end;
sided_pval = ones(1,nv0);
sided_pval(maskvtx) = spval;
end


%% AUXILIAR function
function [spval] = sidedPval(pval,tail)
spval = abs(pval);
if tail==-1
    spval(pval<0) = spval(pval<0)*0.5;
    spval(pval>0) = 1-0.5*spval(pval>0);
elseif tail==1
    spval(pval>0) = spval(pval>0)*0.5;
    spval(pval<0) = 1-0.5*spval(pval<0);
elseif tail~=0
    error('Tail must be -1 for left-sided or 0 for two-sided or 1 for right-sided hypothesis testing');
end
end
