function Xptm_run = fast_polytrendmtx(run,ntrs,nruns,order)
% Xptm = fast_polytrendmtx(run,ntrs,nruns,order)
%

if(nargin ~= 4)
  msg = 'USAGE: Xptm = fast_polytrendmtx(run,ntrs,nruns,order)';
  qoe(msg);error(msg);
end

Xptm = ones(ntrs,1);
t = [0:ntrs-1]'; %'
for n = 1:order
  r0 = t.^n;
  M = eye(ntrs) - Xptm*inv(Xptm'*Xptm)*Xptm';
  r = M*r0;
  r = r/std(r);
  Xptm = [Xptm r];
end

Xptm_run = zeros(ntrs,nruns*(order+1));
n1 = (run-1)*(order+1) + 1;
n2 = n1 + order;
Xptm_run(:,n1:n2) = Xptm;

return;
