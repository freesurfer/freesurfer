function p = fast_z2p(z)
% p = fast_z2p(z)
% converts z values into p values (ie significances)
%

p = erfc(abs(z(:))/sqrt(2))/2;
p = reshape(p,size(z)).*sign(z);
ind = find(z==0);
p(ind) = 1;

return;




