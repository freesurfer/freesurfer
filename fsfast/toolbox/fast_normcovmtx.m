function NormCovMtx = fast_normcovmtx(CovMtx)

d = diag(CovMtx);
ind = find(d==0);
d(ind) = 10^10;
NormCovMtx = CovMtx./sqrt(d*d'); %'

return;
