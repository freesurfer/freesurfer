function dim = fast_estdimsvd(s,pvsthresh)
% dim = fast_estdimsvd(s,<pvsthresh>)

if(nargin ~= 1 & nargin ~= 2)
  fprintf('dim = fast_estdimsvd(s,<pvsthresh>)\n');
  return;
end

if(nargin ~= 2) pvsthresh = 0; end

nf = size(s,1);
ds = diag(s);
pvs = 100*ds/sum(ds);

y = randn(nf,10*nf);
My = y*y'; %'
[uy sy blah] = svd(My);
dsy = diag(sy);
pvsy = 100*dsy/sum(dsy);

d = 100*(pvs-pvsy)./pvsy;

dim = max(find(d > pvsthresh));

%nn = 2:20;
%plot(nn,pvsy(nn),nn,pvs(nn));
%keyboard

return;

















