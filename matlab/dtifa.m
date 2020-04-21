function fa = dtifa(T)
% fa = dtifa(T)
% Compute fractional anisotropy of a 3x3 tensor

if(nargin ~= 1)
  fprintf('fa = dtifa(T)\n');
  return;
end

e = eig(T);
emn = mean(e);
d = e - emn;
sse = sum(d.^2);
ss = sum(e.^2);

fa = sqrt(1.5*sse/ss);

return;
