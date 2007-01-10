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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
%    $Revision: 1.3 $
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
