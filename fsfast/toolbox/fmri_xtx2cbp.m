function CBP = xtx2cbp(XtX,nNNCond,nTP)
%
% CBP = xtx2cbp(XtX,nNNCond)
%
%
%
%


%
% fmri_xtx2cbp.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:07 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

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


