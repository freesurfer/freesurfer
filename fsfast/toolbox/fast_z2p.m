function p = fast_z2p(z)
% p = fast_z2p(z)
% converts z values into p values (ie significances)
%

zsize = size(z);
z = reshape1d(z);

p = erfc(z/sqrt(2))/2;
p = reshape(p,zsize);

return;




