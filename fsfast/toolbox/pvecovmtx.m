function [pve, s] = pvecovmtx(m)
% [pve s] = pvecovmtx(m)
%
% Percent Variance Explained.

if(nargin ~= 1) 
  msg = 'USAGE: [pve s] = pvecovmtx(m)';
  qoe(msg); error(msg);
end

s = svd(m);

pve = 100*cumsum(s)/sum(s);

return;
