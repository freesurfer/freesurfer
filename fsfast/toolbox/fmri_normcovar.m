function Cnorm = normcovar(C)
%
% Cnorm = normcovar(C)
%
%
%
%


%
% fmri_normcovar.m
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

nC    = size(C,1);
nRuns = size(C,3);

Cnorm = ones(size(C));

for r = 1: nRuns,
  c = C(:,:,r);
  d = diag(c);
  utri = triu(c,1);
  ltri = tril(c,-1);
  tri  = (utri + ltri')/2; %'
  dm = repmat(d, [1 nC]);
  ntri = tri ./ dm;
  Cnorm(:,:,r) = ntri + ntri' + eye(nC); %'
end


return;
