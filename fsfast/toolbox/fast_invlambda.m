function invLambda = fast_invlambda(ErrCovMtx,nEigen)
% invLambda = fast_invlambda(ErrCovMtx,nEigen)

if(nargin ~= 2)
  msg = 'USAGE: invLambda = fast_invlambda(ErrCovMtx,nEigen)';
  qoe(msg);error(msg);
end

NormErrCovMtx = fast_normcovmtx(ErrCovMtx);

[u s v] = svd(NormErrCovMtx);

nn = [1:nEigen];
ds = diag(s);
ds2 = ds(nn);
ids = 1./ds2;
is = diag(ids);
invLambda = u(:,nn) * is * u(:,nn)'; %'

return;
