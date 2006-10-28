function p = fast_z2p(z)
% p = fast_z2p(z)
% converts z values into p values (ie, the area under the curve 
% to the right of the z value). This is a one-tailed test (ie,
% the p-value does not take the sign of the z).
%
% $Id: fast_z2p.m,v 1.3 2006/10/28 18:56:11 greve Exp $

p = erfc((z(:))/sqrt(2))/2;
p = reshape(p,size(z));

% Do this if you are signing the z
%p = p.*sign(z);
%ind = find(z==0);
%p(ind) = 1;

return;




