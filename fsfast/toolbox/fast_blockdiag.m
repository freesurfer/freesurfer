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
% $Id: fast_blockdiag.m,v 1.2 2004/04/12 01:03:47 greve Exp $


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
