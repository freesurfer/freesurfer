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
%    $Date: 2011/03/02 00:04:03 $
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
