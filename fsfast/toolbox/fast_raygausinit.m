function params = fast_raygausinit(y)
% params = fast_raygausinit(y)
% params = [alpha Rmu Gmu Gstd];
%
% $Id: fast_raygausinit.m,v 1.1 2006/04/25 01:31:40 greve Exp $

if(nargin ~= 1)
  fprintf('params = fast_raygausinit(y)\n');
  return;
end

md = mode(y);
indray = find(y < md);
indgaus = find(y >= md);
alpha = .5;
Rmu = mean(y(indray));
Gmu = mean(y(indgaus));
Gstd = std(y(indgaus));

params = [alpha Rmu Gmu Gstd];

return;
