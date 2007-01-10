function B = fast_blockident(M,nr,nc)
% B = fast_blockident(M,nr,nc)
%
% Puts M on diagonal. There will be nr X nc submatrices
% the same size as M, those off the diagonal will be 0.
% This is refered to as the "block idenity matrix" because
% each matrix on the diagonal is idenitcal.


%
% fast_blockident.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
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

[Mnr Mnc] = size(M);
Bnr = nr*Mnr;
Bnc = nc*Mnc;

B0 = zeros(Mnr,Bnc);
B0(:,1:Mnc) = M;

B = B0;
for n = 2:nr;
  Btmp = fast_mshift(B0, [0 (n-1)*Mnc], 0);
  B = [B; Btmp];
end

return;
