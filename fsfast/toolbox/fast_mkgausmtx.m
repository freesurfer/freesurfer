function G = fast_mkgausmtx(nstddev,len)
% G = fast_mkgausmtx(nstddev,len)
% Create a len-by-len guassian filter matrix 

G = [];

if(nargin ~= 2)
  fprintf('G = fast_mkgausmtx(nstddev,len)\n');
  return;
end

for n = 1:len
  g = fast_gaussian(n,nstddev,len);
  % g = g/sqrt(sum(g.^2));
  g = g/sum(abs(g));
  G = [G; g];
end

return;
