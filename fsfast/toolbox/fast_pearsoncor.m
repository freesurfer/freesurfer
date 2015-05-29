function rho = fast_pearsoncor(x,y)
% rho = fast_pearsoncor(x,y)

[nt nv] = size(x);

dx = x - repmat(mean(x),[nt 1]);
dy = y - repmat(mean(y),[nt 1]);

rho = sum(dx.*dy)./sqrt(sum(dx.*dx).*sum(dy.*dy));

return;

