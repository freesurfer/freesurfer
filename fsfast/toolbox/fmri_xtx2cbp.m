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
%    $Date: 2007/01/10 22:02:34 $
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


