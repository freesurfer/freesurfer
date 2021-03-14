function NormCovMtx = fast_normcovmtx(CovMtx)


%
% fast_normcovmtx.m
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

d = diag(CovMtx);
ind = find(d==0);
d(ind) = 10^10;
NormCovMtx = CovMtx./sqrt(d*d'); %'

return;
