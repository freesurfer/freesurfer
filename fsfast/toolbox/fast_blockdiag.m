function B = fast_blockdiag(M)
% B = fast_blockdiag(M)
%
% M is nr X nc X nd. Each submatrix is placed on the
% the block diagonal such that the final matrix B will
% consisist of nd^2 submatrices, each nr X nc. The 
% off-diagonal blocks will be zero.
%
% See also fast_blockdiag2
%
%


%
% fast_blockdiag.m
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


B = [];

if(nargin ~= 1)
  fprintf('USAGE: B = fast_blockdiag(M)\n');
  return;
end

[Mnr Mnc Mnd] = size(M);
Bnr = Mnd*Mnr;
Bnc = Mnd*Mnc;

B = zeros(Bnr,Bnc);

for n = 1:Mnd
  rmin = (n-1)*Mnr + 1;
  rmax = rmin + Mnr -1;
  rind = rmin:rmax;

  cmin = (n-1)*Mnc + 1;
  cmax = cmin + Mnc -1;
  cind = cmin:cmax;

  B(rind,cind) = M(:,:,n);
end


return;
