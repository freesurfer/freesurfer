function u = fmri_norm(v,order)
% u = fmri_norm(v, <order>)
% 
% normalizes the columns of v
%
% '$Id: fmri_norm.m,v 1.1 2003/03/04 20:47:40 greve Exp $'

if(nargin == 0)
  msg = 'USAGE: u = fmri_norm(v, <order>)';
  qoe(msg);error(msg);
end

if(size(v,1)==1)
  u = ones(size(v));
  return;
end

if(nargin == 1 ) order = 1 ; end

f = (sum(v.^order)).^(1/order);
ind = find(f==0);
f(ind) = 10^10;
u = v ./ repmat(f,[size(v,1) 1]);

return;
