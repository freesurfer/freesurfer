function CBP = xtx2cbp(XtX,nNNCond,nTP)
%
% CBP = xtx2cbp(XtX,nNNCond)
%
%
%
% $Id: fmri_xtx2cbp.m,v 1.1 2003/03/04 20:47:40 greve Exp $

CBP  = zeros(size(XtX));
[nXtX dummy nRuns] = size(XtX);
nH = nXtX/nNNCond;

d = diag(XtX);
nPerCond = d(1:nH:nXtX);

npcm = zeros(nXtX);

for c1 = 1 : nNNCond,
  for c2 = 1 : nNNCond,

     r =  [ 1 : nH ]' + (c1 - 1) * nH;
     c =  [ 1 : nH ]  + (c2 - 1) * nH;

     rr = repmat(r, [1 nH]);
     cc = repmat(c, [nH 1]);
     i = (cc-1)*nXtX + rr;

     if(c1 == c2) n = nPerCond(c1);
     else         n = nTP;
     end
     
     CBP(i) = XtX(i)/n;

  end
end

return;


