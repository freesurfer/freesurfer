function [cnd, mineig, S] = fast_acfcond(acf,taper)
% [cnd mineig S] = fast_acfcond(acf,<taper>)

cnd = [];
mineig = [];
S = [];

if(nargin ~= 1 & nargin ~= 2)
  fprintf('[cnd mineig S] = fast_acfcond(acf,<taper>)\n');
  return;
end

if(exist('taper') ~= 1) taper = []; end
if(~isempty(taper))
  acf = acf .* taper;
end

S = toeplitz(acf);
mineig = min(eig(S));
cnd = cond(S);

return;







