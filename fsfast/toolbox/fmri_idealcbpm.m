function [CBP, XtX] = IdealCBP(nPerCond, CBOrder)
%
% [CBP XtX]= IdealCBP(nPerCond, CBOrder)
%
% Computes Ideal Counterbalancing Probability Matrix
% Assumes fixation is condition 0.
%
% $Id: fmri_idealcbpm.m,v 1.1 2003/03/04 20:47:39 greve Exp $


nCond   = length(nPerCond); % number of conditions
nNNCond = nCond - 1;        % number of non-null conditions
nStim   = sum(nPerCond);    % total number of Stim for all conditions
PrCond  = nPerCond(2:nCond)/nStim;   % probability of stim from non-null condition

nCBP = nNNCond * CBOrder;
CBP  = zeros(nCBP);
XtX  = zeros(nCBP);

for c1 = 1 : nNNCond,
  for c2 = 1 : nNNCond,

     r =  [ 1 : CBOrder ]' + (c1 - 1) * CBOrder;
     c =  [ 1 : CBOrder ]  + (c2 - 1) * CBOrder;

     rr = repmat(r, [1 CBOrder]);
     cc = repmat(c, [CBOrder 1]);
     i = (cc-1)*nCBP + rr;

     if(c1 == c2)
        m = PrCond(c1) * ones(CBOrder);
        n = nPerCond(c1+1);
        d = eye(CBOrder); 
     else
        m  = PrCond(c1) * PrCond(c2) * ones(CBOrder);
        % n = nPerCond(c1+1) +  nPerCond(c2+1);
        n = nStim;
        d =  zeros(CBOrder);
     end
     
     CBP(i) = d + triu(m,1) + tril(m,-1);
     XtX(i) = n*CBP(i);

  end
end

return;


