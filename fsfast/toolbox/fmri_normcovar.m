function Cnorm = normcovar(C)
%
% Cnorm = normcovar(C)
%
%
% $Id: fmri_normcovar.m,v 1.1 2003/03/04 20:47:40 greve Exp $
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
